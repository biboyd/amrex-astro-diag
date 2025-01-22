//
// This program takes a 3-d cartesian plotfile and calculates
// the A=23 beta decay and electron capture rates as well as related nu loss rates
// and thermal nu loss rates
//
//
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


#include <table_rates.H>
#include <burn_type.H>
#include <sneut5.H>
//#include <AMReX_MultiFabUtil.H>
//#include <AMReX_MultiFab.H>


using namespace amrex;
using namespace rate_tables;

void get_nu_losses() {
    // get name for input and output plotfiles
    std::string pltfile(diag_rp::plotfile);

    if (pltfile.empty()) {
        std::cout << "no plotfile specified" << std::endl;
        std::cout << "use: diag.plotfile=plt00000 (for example)" << std::endl;
        amrex::Error("no plotfile");
    }

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    std::string outfile = "nu_loss." + std::filesystem::path(pltfile).filename().string();

    //get plotdata   
    PlotFileData pf(pltfile);
    const int nlevs = pf.finestLevel() + 1;

    // init vectors needed to writeout
    Vector<std::string> varnames = pf.varNames();
    Vector<int> lev_steps(nlevs);
    Vector<IntVect> ref_ratio(nlevs);
    Vector<MultiFab> new_mfs(nlevs);
    Vector<Geometry> geoms(nlevs);


    // find variable indices -- we want density, temperature, and species.
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    // the plotfile can store either (rho X) or just X alone.  Here we'll assume
    // that we have just X alone

    int rho_comp = get_dens_index(varnames);
    int temp_comp = get_temp_index(varnames);
    int pres_comp = get_pres_index(varnames);
    int spec_comp = get_spec_index(varnames);
    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> out_varnames;
    out_varnames.push_back("rho");
    out_varnames.push_back("X(Na23)");
    out_varnames.push_back("X(Ne23)");
    out_varnames.push_back("A23_electron_capture_rate");
    out_varnames.push_back("A23_beta_decay_rate");   
    out_varnames.push_back("A23_electron_capture_nu_loss");
    out_varnames.push_back("A23_beta_decay_nu_loss");
    out_varnames.push_back("thermal_nu_loss");
       
    // init the rhs. reaction stuff
    //init_tabular();

    for (int ilev = pf.finestLevel(); ilev >= 0; --ilev) {
        // read plotfile data
        lev_steps[ilev] = pf.levelStep(ilev);
        ref_ratio[ilev] = IntVect(pf.refRatio(ilev), 
                                  pf.refRatio(ilev), 
                                  pf.refRatio(ilev));

        const MultiFab mf = pf.get(ilev);
        const MultiFab temp_mf = pf.get(ilev, varnames[temp_comp]);
        const MultiFab press_mf = pf.get(ilev, varnames[pres_comp]);
        const MultiFab rho_mf = pf.get(ilev, varnames[rho_comp]);
        //const MultiFab spec_mf = pf.get(ilev, varnames[spec_comp]);

        // define new mf with 5 componenets
        // 2 rates, 2 nu energy rate, 1 thermal energy rate
        new_mfs[ilev].define(mf.boxArray(),
                        mf.DistributionMap(),
                        static_cast<int>(out_varnames.size()),
                        mf.nGrow());

        // construct geometry from plotfile
        const amrex::RealBox rb(pf.probLo(), pf.probHi());
        const Array<int,3> is_per({0, 0, 0});
        geoms[ilev].define(pf.probDomain(ilev), rb, pf.coordSys(), is_per);

        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi){
            // find bounds
            const auto tileBox = mfi.tilebox();
          
            // load arrays
            Array4<const Real> const& temp_arr = temp_mf.const_array(mfi);
            Array4<const Real> const& pres_arr = press_mf.const_array(mfi);
            Array4<const Real> const& rho_arr = rho_mf.const_array(mfi);
            Array4<const Real> const& X_arr = mf.const_array(mfi, spec_comp);
            Array4<Real> const& new_arr = new_mfs[ilev].array(mfi);

            // loop over all cells
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                // initialize EOS
                eos_t eos_state;
                burn_t burn_state;

                eos_state.T    = temp_arr(i, j, k);
                eos_state.p    = pres_arr(i, j, k);
                eos_state.rho  = rho_arr(i, j, k);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = X_arr(i, j, k, comp);
                }
                eos(eos_input_rt, eos_state); // use rt instead?
     
                int ine23 = network_spec_index("neon-23");
                int ina23 = network_spec_index("sodium-23");
             
                composition(eos_state);
                eos_to_burn(eos_state, burn_state);
             
                amrex::Real rate, drate_dt, edot_nu, edot_gamma;
                amrex::Real rhoy = burn_state.rho * burn_state.y_e;


#ifndef SKIP_ECAP
                tabular_evaluate(j_Na23_Ne23_meta, j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);
                Real r_ecap = rate;
                Real specific_energy_ecap = C::Legacy::n_A * burn_state.xn[ina23]/23 * (edot_nu + edot_gamma);
#else
                Real r_ecap = 0.;
                Real specific_energy_ecap = 0.;
#endif

#ifndef SKIP_BETA 
                tabular_evaluate(j_Ne23_Na23_meta, j_Ne23_Na23_rhoy, j_Ne23_Na23_temp, j_Ne23_Na23_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);
                Real r_beta = rate;
                Real specific_energy_beta = C::Legacy::n_A * burn_state.xn[ine23]/23 * (edot_nu + edot_gamma);
#else
                Real r_beta = 0.;
                Real specific_energy_beta = 0.;
#endif
                Real xr_ecap = burn_state.xn[ina23] * r_ecap;
                Real xr_beta = burn_state.xn[ine23] * r_beta;


                //thermal neutrino loss   
                Real sneut, dsneutdt, dsneutdd, dsnuda, dsnudz;

                constexpr int do_T_derivatives = 0;
                sneut5<do_T_derivatives>(burn_state.T, burn_state.rho, burn_state.abar, burn_state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

                //save values
                new_arr(i, j, k, 0) = rho_arr(i, j, k);
                new_arr(i, j, k, 1) = X_arr(i, j, k, ina23); 
                new_arr(i, j, k, 2) = X_arr(i, j, k, ine23);
                new_arr(i, j, k, 3) = xr_ecap;
                new_arr(i, j, k, 4) = xr_beta;
                new_arr(i, j, k, 5) = specific_energy_ecap;
                new_arr(i, j, k, 6) = specific_energy_beta;
                new_arr(i, j, k, 7) = sneut;
            });
        }
    }
                   
    //const Vector<const MultiFab* > write_mfs = new_mfs;
    const Vector<Geometry> write_geoms = geoms;

    // write out with gradient fields
    WriteMultiLevelPlotfile(outfile,
                            nlevs,
                            GetVecOfConstPtrs(new_mfs),
                            out_varnames,
                            write_geoms,
                            pf.time(),
                            lev_steps,
                            ref_ratio);

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

    get_nu_losses();
    amrex::Finalize();
}