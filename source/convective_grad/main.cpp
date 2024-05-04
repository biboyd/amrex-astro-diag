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

    std::string outfile = "convgrad." +
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
    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;
    gvarnames.push_back("del");
    gvarnames.push_back("del_ad");
    gvarnames.push_back("del_ledoux");

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

    // get center if spherical

    Array<Real, AMREX_SPACEDIM> center;
    auto const probLo = pf.probLo();
    auto const probHi = pf.probHi();

    if (diag_rp::spherical){
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            center[idim] = 0.5_rt * (probHi[idim] - probLo[idim]);
        }
    }

    // we need both T and P constructed with ghost cells

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), static_cast<int>(gvarnames.size()), 0);

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

            // temperature and pressure with ghost cells
            auto const& T = temp_mf.const_array(mfi);
            auto const& P = pres_mf.const_array(mfi);
            auto const& X = species_mf.const_array(mfi);

            // all of the data without ghost cells
            const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                // calc position if spherical
                Real xpos, ypos, zpos;
                if (diag_rp::spherical){
                    xpos = probLo[0] + dx[0] * (Real(i) + 0.5_rt) - center[0];
                    ypos = probLo[1] + dx[1] * (Real(j) + 0.5_rt) - center[1];
                    zpos = probLo[2] + dx[2] * (Real(k) + 0.5_rt) - center[2];
                }

                // first dlog T / dlog P actual -- we assume that the last
                // dimension is the vertical (plane-parallel)

                if (!diag_rp::spherical) {

                    if (ndims == 1) {
                        // x is the vertical
                        Real dp = P(i+1,j,k) - P(i-1,j,k);
                        if (dp != 0.0) {
                            ga(i,j,k,0) = (T(i+1,j,k) - T(i-1,j,k)) / dp * (P(i,j,k) / T(i,j,k));
                        } else {
                            ga(i,j,k,0) = 0.0;
                        }

                    } else if (ndims == 2) {
                        // y is the vertical
                        Real dp = P(i,j+1,k) - P(i,j-1,k);
                        if (dp != 0.0) {
                            ga(i,j,k,0) = (T(i,j+1,k) - T(i,j-1,k)) / dp * (P(i,j,k) / T(i,j,k));
                        } else {
                            ga(i,j,k,0) = 0.0;
                        }

                    } else {
                        // z is the vertical
                        Real dp = P(i,j,k+1) - P(i,j,k-1);
                        if (dp != 0.0) {
                            ga(i,j,k,0) = (T(i,j,k+1) - T(i,j,k-1)) / dp * (P(i,j,k) / T(i,j,k));
                        } else {
                            ga(i,j,k,0) = 0.0;
                        }
                    }
                } else {
                    //spherical case
                    if (ndims == 1) {
                        // x is the vertical
                        // same as plane parallel
                        Real dp = P(i+1,j,k) - P(i-1,j,k);
                        if (dp != 0.0) {
                            ga(i,j,k,0) = (T(i+1,j,k) - T(i-1,j,k)) / dp * (P(i,j,k) / T(i,j,k));
                        } else {
                            ga(i,j,k,0) = 0.0;
                        }

                    } else if (ndims == 2) {
                        // r is from x and y
                        Real dp = (xpos / dx[0]) * (P(i+1,j,k) - P(i-1,j,k))
                                + (ypos / dx[1]) * (P(i,j+1,k) - P(i,j-1,k));

                        if (dp != 0.0) {
                            Real dT = (xpos / dx[0]) * (T(i+1,j,k) - T(i-1,j,k))
                                    + (ypos / dx[1]) * (T(i,j+1,k) - T(i,j-1,k));

                            ga(i,j,k,0) = (dT / dp) * (P(i,j,k) / T(i,j,k));
                        } else {
                            ga(i,j,k,0) = 0.0;
                        }

                    } else {
                       // r is from x, y, and z
                        Real dp = (xpos / dx[0]) * (P(i+1,j,k) - P(i-1,j,k))
                                + (ypos / dx[1]) * (P(i,j+1,k) - P(i,j-1,k))
                                + (zpos / dx[2]) * (P(i,j,k+1) - P(i,j,k-1));

                        if (dp != 0.0) {
                            Real dT = (xpos / dx[0]) * (T(i+1,j,k) - T(i-1,j,k))
                                    + (ypos / dx[1]) * (T(i,j+1,k) - T(i,j-1,k))
                                    + (zpos / dx[2]) * (T(i,j,k+1) - T(i,j,k-1));

                            ga(i,j,k,0) = (dT / dp) * (P(i,j,k) / T(i,j,k));

                        } else {
                            ga(i,j,k,0) = 0.0;
                        }
                    }
                }

                // now del_ad.  We'll follow HKT Eq. 3.96, 3.97

                eos_t eos_state;

                eos_state.rho = fab(i,j,k,dens_comp);
                eos_state.T = fab(i,j,k,temp_comp);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j,k,n);
                }
                eos(eos_input_rt, eos_state);

                Real chi_T = eos_state.dpdT * eos_state.T / eos_state.p;

                ga(i,j,k,1) = eos_state.p * chi_T / (eos_state.gam1 * eos_state.rho * eos_state.T * eos_state.cv);


                // del_ledoux = del_ad + B, where B is the composition term
                // We calculate it like MESA, Paxton+ 2013 Equation 8
                // but we do a centered difference

                Real lnP_plus{0.0};  // pressure "above"
                Real lnP_minus{0.0};  // pressure "below"

                Real lnPalt_plus{0.0};  // pressure with "above" species
                Real lnPalt_minus{0.0};  // pressure with "below" species

                if (! diag_rp::spherical){

                    if (ndims == 1) {
                        // x is the vertical

                        lnP_plus = std::log(P(i+1,j,k));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i+1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus = std::log(eos_state.p);

                        lnP_minus = std::log(P(i-1,j,k));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i-1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus = std::log(eos_state.p);

                    } else if (ndims ==2 ) {
                        // y is the vertical

                        lnP_plus = std::log(P(i,j+1,k));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j+1,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus = std::log(eos_state.p);

                        lnP_minus = std::log(P(i,j-1,k));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j-1,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus = std::log(eos_state.p);

                    } else {
                        // z is the vertical

                        lnP_plus = std::log(P(i,j,k+1));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j,k+1,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus = std::log(eos_state.p);

                        lnP_minus = std::log(P(i,j,k-1));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j,k-1,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus = std::log(eos_state.p);
                    }
                } else{
                    //spherical case

                    if (ndims == 1) {
                        // x is the vertical
                        // same as plane-parallel

                        lnP_plus = std::log(P(i+1,j,k));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i+1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus = std::log(eos_state.p);

                        lnP_minus = std::log(P(i-1,j,k));
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i-1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus = std::log(eos_state.p);

                    } else if (ndims ==2 ) {
                        // r is made of x and y

                        // actual
                        lnP_plus = (xpos / dx[0]) * std::log(P(i+1,j,k))
                                 + (ypos / dx[1]) * std::log(P(i,j+1,k));

                        lnP_minus = (xpos / dx[0]) * std::log(P(i-1,j,k))
                                  + (ypos / dx[1]) * std::log(P(i,j-1,k));

                        //alternate
                        //plus - x
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i+1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus += (xpos / dx[0]) * std::log(eos_state.p);

                        //plus - y
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j+1,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus += (ypos / dx[1]) * std::log(eos_state.p);

                        //minus - x
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i-1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus += (xpos / dx[0]) * std::log(eos_state.p);

                        //minus - y
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j-1,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus += (ypos / dx[1]) * std::log(eos_state.p);

                    } else {
                        // r is made of x, y and z

                        // actual
                        lnP_plus = (xpos / dx[0]) * std::log(P(i+1,j,k))
                                 + (ypos / dx[1]) * std::log(P(i,j+1,k))
                                 + (zpos / dx[2]) * std::log(P(i,j,k+1));

                        lnP_minus = (xpos / dx[0]) * std::log(P(i-1,j,k))
                                  + (ypos / dx[1]) * std::log(P(i,j-1,k))
                                  + (zpos / dx[2]) * std::log(P(i,j,k-1));

                        //alternate
                        //plus - x
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i+1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus += (xpos / dx[0]) * std::log(eos_state.p);

                        //plus - y
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j+1,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus += (ypos / dx[1]) * std::log(eos_state.p);

                        //plus - z
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j,k+1,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_plus += (zpos / dx[2]) * std::log(eos_state.p);

                        //minus - x
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i-1,j,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus += (xpos / dx[0]) * std::log(eos_state.p);

                        //minus - y
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j-1,k,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus += (ypos / dx[1]) * std::log(eos_state.p);

                        //minus - z
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = X(i,j,k-1,n);
                        }
                        eos(eos_input_rt, eos_state);
                        lnPalt_minus += (zpos / dx[2]) * std::log(eos_state.p);
                    }
                }

                // chi_T still has the old (correct) value for i,j,k

                Real denom = lnP_plus - lnP_minus;
                Real B{0.0};
                if (denom != 0.0) {
                    B = -1 / chi_T * (lnPalt_plus - lnPalt_minus) / denom;
                }
                ga(i,j,k,2) = ga(i,j,k,1) + B;

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
