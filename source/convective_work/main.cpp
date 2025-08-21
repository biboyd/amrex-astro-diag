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

inline Real get_mu_e(eos_t& state){
    // return electron chem potential from degeneracy param eta
    return state.eta * C::k_B * state.T + C::m_e * C::c_light * C::c_light;
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

    std::string outfile = "convwork." +
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
    int spec_comp = get_spec_index(var_names_pf);
    int radvel_comp = get_radvel_index(var_names_pf); //account for no radvel and calc on the fly
    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;
    gvarnames.push_back("chem_e");
    gvarnames.push_back("d_chem_e");
    gvarnames.push_back("flux_e");
    gvarnames.push_back("eps_conv");

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

    // get center

    Array<Real, AMREX_SPACEDIM> center;
    auto const probLo = pf.probLo();
    auto const probHi = pf.probHi();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
        center[idim] = 0.5_rt * (probHi[idim] - probLo[idim]);
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

        // fill the density and temperature mfs with ghost cells
        // we also need all of the species

        MultiFab temp_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab dens_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab species_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), NumSpec, ng);

        if (ilev == 0) {

            // temperature
            {
                MultiFab smf = pf.get(ilev, var_names_pf[temp_comp]);
                FillPatchSingleLevel(temp_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // density 
            {
                MultiFab smf = pf.get(ilev, var_names_pf[dens_comp]);
                FillPatchSingleLevel(dens_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
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

            // density 
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[dens_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[dens_comp]);
                FillPatchTwoLevels(dens_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
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

            // temperature and density with ghost cells
            auto const& rho = dens_mf.const_array(mfi);
            auto const& T = temp_mf.const_array(mfi);
            auto const& X = species_mf.const_array(mfi);

            // all of the data without ghost cells
            const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                // calc position 
                Real xpos, ypos, zpos, rpos;
                xpos = probLo[0] + dx[0] * (Real(i) + 0.5_rt) - center[0];
                ypos = probLo[1] + dx[1] * (Real(j) + 0.5_rt) - center[1];
                zpos = probLo[2] + dx[2] * (Real(k) + 0.5_rt) - center[2];
                rpos = std::sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

                // calc Flux of electrons
                eos_t eos_state;

                eos_state.rho = rho(i,j,k);
                eos_state.T = T(i,j,k);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j,k,n);
                }
                eos(eos_input_rt, eos_state);

                //write flux and mu_e out
                ga(i, j, k, 0) = get_mu_e(eos_state);

                Real F_e = C::n_A * eos_state.rho * eos_state.y_e * fab(i,j,k,radvel_comp);
                ga(i, j, k, 2) = F_e;


                // Now dlog mu_e / dr 
                // Note: don't need m_e c^2 term b/c will cancel out
                //dmu_dx
                Real dmu_dx{0.}, dmu_dy{0.}, dmu_dz{0.}, dmu_dr{0.}; 

                //x
                // add i+1
                eos_state.rho = rho(i+1,j,k);
                eos_state.T = T(i+1,j,k);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i+1,j,k,n);
                }
                eos(eos_input_rt, eos_state);
                dmu_dx += eos_state.eta * eos_state.T * C::k_B;

                // sub i-1
                eos_state.rho = rho(i-1,j,k);
                eos_state.T = T(i-1,j,k);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i-1,j,k,n);
                }
                eos(eos_input_rt, eos_state);
                dmu_dx -= eos_state.eta * eos_state.T * C::k_B;

                // div 2x
                dmu_dx /= 2 * dx[0];

                //y
                // add j+1
                eos_state.rho = rho(i,j+1,k);
                eos_state.T = T(i,j+1,k);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j+1,k,n);
                }
                eos(eos_input_rt, eos_state);
                dmu_dy += eos_state.eta * eos_state.T * C::k_B;

                // sub j-1
                eos_state.rho = rho(i,j-1,k);
                eos_state.T = T(i,j-1,k);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j-1,k,n);
                }
                eos(eos_input_rt, eos_state);
                dmu_dy -= eos_state.eta * eos_state.T * C::k_B;

                // div 2y
                dmu_dy /= 2 * dx[1];

                //z
                //add k+1
                eos_state.rho = rho(i,j,k+1);
                eos_state.T = T(i,j,k+1);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j,k+1,n);
                }
                eos(eos_input_rt, eos_state);
                dmu_dz += eos_state.eta * eos_state.T * C::k_B;

                //sub k-1
                eos_state.rho = rho(i,j,k-1);
                eos_state.T = T(i,j,k-1);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j,k-1,n);
                }
                eos(eos_input_rt, eos_state);
                dmu_dz -= eos_state.eta * eos_state.T * C::k_B;

                // div 2z
                dmu_dz /= 2 * dx[2];

                // r and save
                dmu_dr = (xpos*dmu_dx + ypos*dmu_dy + zpos*dmu_dz)/rpos;
                ga(i,j,k,1) = dmu_dr;
                ga(i,j,k,3) = F_e * dmu_dr;
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
