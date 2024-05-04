//
// Process a sedov problem to produce rho, u, and p as a
// function of r, for comparison to the analytic solution.
//
#include <iostream>
// #include <stringstream>
#include <regex>
#include <string>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <amrex_astro_util.H>

#include <extern_parameters.H>

#include <network.H>
#include <eos.H>

using namespace amrex;


int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv);

    // initialize the runtime parameters

    init_extern_parameters();

    // initialize C++ Microphysics

    eos_init(diag_rp::small_temp, diag_rp::small_dens);
    network_init();

    // timer for profiling

    BL_PROFILE_VAR("main()", pmain);

    // read the plotfile metadata

    PlotFileData pf(diag_rp::plotfile);

    int fine_level = pf.finestLevel();
    const int dim = pf.spaceDim();

    // find variable indices -- we want density, temperature, and species.
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    // the plotfile can store either (rho X) or just X alone.  Here we'll assume
    // that we have just X alone

    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = get_dens_index(var_names_pf);
    int temp_comp = get_temp_index(var_names_pf);
    int spec_comp = get_spec_index(var_names_pf);

    // we will use a mask that tells us if a zone on the current level
    // is covered by data on a finer level.

    for (int ilev = 0; ilev <= fine_level; ++ilev) {

        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                pf.boxArray(ilev+1), ratio);

            const MultiFab& lev_data_mf = pf.get(ilev);

            for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.ok()) {
                    const auto& m = mask.array(mfi);
                    const auto& fab = lev_data_mf.array(mfi);
                    const auto lo = amrex::lbound(bx);
                    const auto hi = amrex::ubound(bx);

                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {

                                // if we only wanted to work on the zones not covered
                                // by finer grids (e.g., for computing a global sum)
                                // we'd uncomment this mask check

                                //if (m(i,j,k) == 0) { // not covered by fine

                                // compute the coordinate of the current zone

                                eos_t eos_state;

                                eos_state.rho = fab(i,j,k,dens_comp);
                                eos_state.T = fab(i,j,k,temp_comp);
                                for (int n = 0; n < NumSpec; ++n) {
                                    eos_state.xn[n] = fab(i,j,k,spec_comp+n);
                                }

                                eos(eos_input_rt, eos_state);

                                // } // mask

                            }
                        }
                    }

                } // bx.ok()

            } // MFIter

        } else {
            // this is the finest level

            const MultiFab& lev_data_mf = pf.get(ilev);

            for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.ok()) {
                    const auto& fab = lev_data_mf.array(mfi);
                    const auto lo = amrex::lbound(bx);
                    const auto hi = amrex::ubound(bx);

                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {

                                eos_t eos_state;

                                eos_state.rho = fab(i,j,k,dens_comp);
                                eos_state.T = fab(i,j,k,temp_comp);
                                for (int n = 0; n < NumSpec; ++n) {
                                    eos_state.xn[n] = fab(i,j,k,spec_comp+n);
                                }

                                eos(eos_input_rt, eos_state);

                            }
                        }
                    }

                } // bx.ok()

            } // MFIter


        }

    } // level loop



    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}

