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

int get_var_idx(const std::string& var_str, const std::vector<std::string>& var_names_pf) {
   //loop over all variables to get idx
   auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), var_str);
   if (idx == var_names_pf.cend()) {
        amrex::Print() << var_str << std::endl;
        amrex::Error("Error: could not find a variable");
        }
   int var_idx = std::distance(var_names_pf.cbegin(), idx);
   return var_idx;
}

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

    std::string outfile = "small." +
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

    Vector<std::string> gvar_names;
    Vector<int> gvar_idxs;
    gvar_names.push_back("velx");
    gvar_names.push_back("vely");
    gvar_names.push_back("velz");
    gvar_names.push_back("vort");
    gvar_names.push_back("rho");
    gvar_names.push_back("h");
    gvar_names.push_back("X(n)");
    gvar_names.push_back("X(H1)");
    gvar_names.push_back("X(He4)");
    gvar_names.push_back("X(C12)");
    gvar_names.push_back("X(C13)");
    gvar_names.push_back("X(C14)");
    gvar_names.push_back("X(N13)");
    gvar_names.push_back("X(N14)");
    gvar_names.push_back("X(O16)");
    gvar_names.push_back("X(O17)");
    gvar_names.push_back("X(O18)");
    gvar_names.push_back("X(F18)");
    gvar_names.push_back("X(F21)");
    gvar_names.push_back("X(Ne20)");
    gvar_names.push_back("X(Ne21)");
    gvar_names.push_back("X(Ne22)");
    gvar_names.push_back("X(Ne23)");
    gvar_names.push_back("X(Na23)");
    gvar_names.push_back("X(Na25)");
    gvar_names.push_back("X(Mg24)");
    gvar_names.push_back("X(Mg25)");
    gvar_names.push_back("omegadot(n)");
    gvar_names.push_back("omegadot(H1)");
    gvar_names.push_back("omegadot(He4)");
    gvar_names.push_back("omegadot(C12)");
    gvar_names.push_back("omegadot(C13)");
    gvar_names.push_back("omegadot(C14)");
    gvar_names.push_back("omegadot(N13)");
    gvar_names.push_back("omegadot(N14)");
    gvar_names.push_back("omegadot(O16)");
    gvar_names.push_back("omegadot(O17)");
    gvar_names.push_back("omegadot(O18)");
    gvar_names.push_back("omegadot(F18)");
    gvar_names.push_back("omegadot(F21)");
    gvar_names.push_back("omegadot(Ne20)");
    gvar_names.push_back("omegadot(Ne21)");
    gvar_names.push_back("omegadot(Ne22)");
    gvar_names.push_back("omegadot(Ne23)");
    gvar_names.push_back("omegadot(Na23)");
    gvar_names.push_back("omegadot(Na25)");
    gvar_names.push_back("omegadot(Mg24)");
    gvar_names.push_back("omegadot(Mg25)");
    gvar_names.push_back("Hnuc");
    gvar_names.push_back("tfromp");
    gvar_names.push_back("p0pluspi");
    gvar_names.push_back("entropy");

    for (auto name_it = gvar_names.begin(); name_it < gvar_names.end(); ++name_it ){
        gvar_idxs.push_back(get_var_idx(*name_it, var_names_pf));
    }

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), static_cast<int>(gvar_names.size()), 0);

        Vector<BCRec> bcr{bcr_default};
    	

        Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                         pf.coordSys(), is_periodic);
        geom.push_back(vargeom);

        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
            (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));


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
    			for (int n=0; n < gvar_idxs.size(); ++n){
                    ga(i, j, k, n) = fab(i, j, k, gvar_idxs[n]);
                }
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

    WriteMultiLevelPlotfile(outfile, nlevs, GetVecOfConstPtrs(gmf), gvar_names,
                            geom, pf.time(), level_steps, ref_ratio);
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv);

    // initialize the runtime parameters

    init_extern_parameters();

    // initialize C++ Microphysics

    main_main();
    amrex::Finalize();
}
