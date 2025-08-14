#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>

#include <extern_parameters.H>

#include <network.H>
#include <eos.H>

#include <amrex_astro_util.H>

using namespace amrex;


void main_main()
{

    std::string pltfile(diag_rp::plotfile);

    if (pltfile.empty()) {
        std::cout << "no plotfile specified" << std::endl;
        std::cout << "use: diag.plotfile=plt00000 (for example)" << std::endl;
        amrex::Error("no plotfile");
    }

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    std::string outfile = "eosderiv." +
        std::filesystem::path(pltfile).filename().string();


    PlotFileData pf(pltfile);

    const int ndims = pf.spaceDim();
    AMREX_ALWAYS_ASSERT(ndims <= AMREX_SPACEDIM);

    const int nlevs = pf.finestLevel() + 1;

    Vector<std::string> varnames;
    varnames = pf.varNames();

    // find variable indices -- we want density, temperature, and species.
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    // the plotfile can store either (rho X) or just X alone.  Here we'll assume
    // that we have just X alone

    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = get_dens_index(var_names_pf);
    int temp_comp = get_temp_index(var_names_pf);
    int pres_comp = get_pres_index(var_names_pf);
    int spec_comp = get_spec_index(var_names_pf);
    int omegadot_comp = get_omegadot_index(var_names_pf);


    // interpret the boundary conditions

    BCRec bcr_default;
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    IntVect ng(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim < ndims) {
            bcr_default.setLo(idim, BCType::hoextrapcc);
            bcr_default.setHi(idim, BCType::hoextrapcc);
        } else {
            bcr_default.setLo(idim, BCType::int_dir);
            bcr_default.setHi(idim, BCType::int_dir);
            is_periodic[idim] = 1;
            ng[idim] = 0;
        }
    }


    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;
    gvarnames.push_back("density");
    gvarnames.push_back("temperature");
    gvarnames.push_back("internal_energy");
    gvarnames.push_back("eta_e");
    gvarnames.push_back("therm_term_h");
    gvarnames.push_back("sigma");

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), static_cast<int>(gvarnames.size()), 0);

        Vector<BCRec> bcr{bcr_default};
    	

        Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                         pf.coordSys(), is_periodic);
        geom.push_back(vargeom);

        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
            (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

        // fill the density and temperature mfs with ghost cells
        // we also need all of the species


        auto const& dx = pf.cellSize(ilev);

        const MultiFab& lev_data_mf = pf.get(ilev);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(lev_data_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();

            // output storage
            auto const& ga = gmf[ilev].array(mfi);


            // all of the data without ghost cells
            const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // calc eos
                eos_extra_t eos_state;

                eos_state.rho = fab(i,j,k,dens_comp);
                eos_state.T = fab(i,j,k,temp_comp);
                eos_state.p = fab(i,j,k,pres_comp);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = fab(i,j,k,spec_comp+n);
                }

                eos(eos_input_rp, eos_state);
            	eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

                //calc sum dhdX. See xhi in S in maestro eqn's
                Real therm_term = 0.;
                for (int n = 0; n < NumSpec; ++n) {
                    therm_term -= eos_xderivs.dhdX[n] * fab(i,j,k,omegadot_comp+n);
                }

				//calc sigma = p_T/(rho cp p_rho)
				Real sigma = eos_state.dpdT / (eos_state.rho*eos_state.cp*eos_state.dpdr);

                //write everything out
                ga(i, j, k, 0) = eos_state.rho;
                ga(i, j, k, 1) = eos_state.T;
                ga(i, j, k, 2) = eos_state.e;
                ga(i, j, k, 3) = eos_state.eta;
                ga(i, j, k, 4) = therm_term;
                ga(i, j, k, 5) = sigma;
            });
        }
    }

    Vector<int> level_steps;
    Vector<IntVect> ref_ratio;
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        level_steps.push_back(pf.levelStep(ilev));
        if (ilev < pf.finestLevel()) {
            ref_ratio.push_back(IntVect(pf.refRatio(ilev)));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ref_ratio[ilev][idim] = 1;
            }
        }
    }

    WriteMultiLevelPlotfile(outfile, nlevs, GetVecOfConstPtrs(gmf), gvarnames,
                            geom, pf.time(), level_steps, ref_ratio);
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv);

    // initialize the runtime parameters

    init_extern_parameters();

    // initialize C++ Microphysics

    eos_init(diag_rp::small_temp, diag_rp::small_dens);
    network_init();

    main_main();
    amrex::Finalize();
}
