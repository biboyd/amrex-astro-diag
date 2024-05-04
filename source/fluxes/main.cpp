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
#include <conductivity.H>

#include <fundamental_constants.H>

#include <amrex_astro_util.H>

using namespace amrex;

inline int
get_vy_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "vely");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find vely component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

inline int
get_vz_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "velz");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find velz component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

inline int
get_dT_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "tpert");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find tpert component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile(diag_rp::plotfile);

    if (pltfile.empty()) {
        std::cout << "no plotfile specified" << std::endl;
        std::cout << "use: diag.plotfile=plt00000 (for example)" << std::endl;
        amrex::Error("no plotfile");
    }

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    std::string outfile = pltfile + "/fluxes";
    std::cout << outfile << std::endl;

    PlotFileData pf(pltfile);

    const int ndims = pf.spaceDim();
    AMREX_ALWAYS_ASSERT(ndims <= AMREX_SPACEDIM);

    const int nlevs = pf.finestLevel() + 1;

    Vector<std::string> varnames;
    varnames = pf.varNames();

    // find variable indices
    // We want:
    // density, temperature, pressure, species
    // vertical velocity, temperature perturbation
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = get_dens_index(var_names_pf);
    int temp_comp = get_temp_index(var_names_pf);
    int pres_comp = get_pres_index(var_names_pf);
    int spec_comp = get_spec_index(var_names_pf);
    int dT_comp = get_dT_index(var_names_pf);

    int v_comp = get_vy_index(var_names_pf);
    if (ndims == 3) {
        // z is the vertical
        int v_comp = get_vz_index(var_names_pf);
    }

    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;
    gvarnames.push_back("Fconv");
    gvarnames.push_back("Fconv_mlt");
    gvarnames.push_back("Fkin");
    gvarnames.push_back("Frad");
    gvarnames.push_back("Fh1");

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

    // we need the variables constructed with ghost cells

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), gvarnames.size(), 0);

        Vector<BCRec> bcr{bcr_default};
        auto is_per = is_periodic;

        Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                         pf.coordSys(), is_per);
        geom.push_back(vargeom);

        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
            (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

        // fill the pressure and temperature mfs with ghost cells
        // we also need all of the species

        MultiFab temp_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab pres_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab species_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), NumSpec, ng);
        MultiFab vy_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab dT_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);

        if (ilev == 0) {

            // temperature
            {
                MultiFab smf = pf.get(ilev, var_names_pf[temp_comp]);
                FillPatchSingleLevel(temp_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // pressure
            {
                MultiFab smf = pf.get(ilev, var_names_pf[pres_comp]);
                FillPatchSingleLevel(pres_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // species
            {
                for (int n = 0; n < NumSpec; ++n) {
                    MultiFab smf = pf.get(ilev, var_names_pf[spec_comp+n]);
                    FillPatchSingleLevel(species_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                         0, n, 1, vargeom, physbcf, 0);
                }
            }

            // vertical velocity
            {
                MultiFab smf = pf.get(ilev, var_names_pf[v_comp]);
                FillPatchSingleLevel(vy_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // dT
            {
                MultiFab smf = pf.get(ilev, var_names_pf[dT_comp]);
                FillPatchSingleLevel(dT_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }



        } else {
            auto* mapper = (Interpolater*)(&cell_cons_interp);

            IntVect ratio(pf.refRatio(ilev-1));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }

            Geometry cgeom(pf.probDomain(ilev-1), RealBox(pf.probLo(),pf.probHi()),
                           pf.coordSys(), is_per);
            PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cphysbcf
                (cgeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

            // temperature
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[temp_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[temp_comp]);
                FillPatchTwoLevels(temp_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // pressure
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[pres_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[pres_comp]);
                FillPatchTwoLevels(pres_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // species
            {
                for (int n = 0; n < NumSpec; ++n) {
                    MultiFab cmf = pf.get(ilev-1, var_names_pf[spec_comp+n]);
                    MultiFab fmf = pf.get(ilev  , var_names_pf[spec_comp+n]);
                    FillPatchTwoLevels(species_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                       {&fmf}, {Real(0.0)}, 0, n, 1, cgeom, vargeom,
                                       cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
                }
            }

            // v
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[v_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[v_comp]);
                FillPatchTwoLevels(vy_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }


            // dT
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[dT_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[dT_comp]);
                FillPatchTwoLevels(dT_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

        }

        auto const& dx = pf.cellSize(ilev);

        const MultiFab& lev_data_mf = pf.get(ilev);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(temp_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();

            // output storage
            auto const& ga = gmf[ilev].array(mfi);

            // temperature with ghost cells
            auto const& T = temp_mf.const_array(mfi);

            // all of the data without ghost cells
            const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                Real dT_dr = 0.0;
                if ( ndims == 2 ) {
                    // y is the vertical
                    dT_dr = (T(i,j+1,k) - T(i,j-1,k)) / (2.0*dx[1]);
                } else {
                    // z is the vertical
                    dT_dr = (T(i,j,k+1) - T(i,j-1,k-1)) / (2.0*dx[2]);
                }


                Real pres = fab(i,j,k,pres_comp);
                Real rho  = fab(i,j,k,dens_comp);
                Real temp = fab(i,j,k,temp_comp);
                Real vel   = fab(i,j,k,v_comp);
                Real delT   = fab(i,j,k,dT_comp);

                // Make EOS
                eos_t eos_state;
                eos_state.rho = rho;
                eos_state.T = temp;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = fab(i,j,k,spec_comp+n);
                }
                eos(eos_input_rt, eos_state);

                conductivity(eos_state);

                // Derive from EOS
                Real cp = eos_state.cp;
                Real Q = temp/rho * eos_state.dpdT/eos_state.dpdr; // dlnd/dlnT = T/d dd/dT = T/d (dP/dT)/(dP/dd) = T/d chi_T/chi_d

                // Other
                // auto grav = GetVarFromJobInfo(pltfile, "maestro.grav_const"); // really slow!!
                // Real g = std::stod(grav);
                // std::cout << g << std::endl;
                //Real Hp = -pres/(rho*g) // g is negative

                // Convective heat flux
                ga(i,j,k,0) = rho * cp * vel * delT;

                // Mixing-length heat flux
                // ga(i,j,k,1) = rho * cp * temp * pow(vel, 3) / (Q * g * Hp);
                //ga(i,j,k,1) = pow(rho,2) * cp * temp * pow(vel,3) / (Q * pres); // doesn't require g
                ga(i,j,k,1) = pow(rho,2) * cp * temp * pow(std::abs(vel), 3) / (Q * pres); // using absolute value of velocity

                // Kinetic flux
                ga(i,j,k,2) = rho * pow(vel,3);

                // Radiative flux
                // conductivity is k = 4*a*c*T^3/(kap*rho)
                // see Microphysics/conductivity/stellar/actual_conductivity.H
                ga(i,j,k,3) = -eos_state.conductivity * dT_dr;

                // Hydrogen flux
                ga(i,j,k,4) = rho * vel * fab(i,j,k,spec_comp+0); // this is rho*v*X, not rho*v*dX

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
