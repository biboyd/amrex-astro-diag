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

    //keep track number of variables for each pairing
    int nA21{0}, nA23{0}, nA25{0};

    Vector<std::string> out_varnames;
    out_varnames.push_back("rho");

//Add A21 nuclei if either is included. 
//Add rates only if react is included
#ifndef SKIP_A21_BETA 
    out_varnames.push_back("X(Ne21)");
    out_varnames.push_back("X(F21)");
    out_varnames.push_back("A21_beta_decay_rate");   
    out_varnames.push_back("A21_beta_decay_nu_loss");
    nA21=4;
#ifndef SKIP_A21_ECAP
    out_varnames.push_back("A21_electron_capture_rate");
    out_varnames.push_back("A21_electron_capture_nu_loss");
    nA21=6;
#endif
#else
#ifndef SKIP_A21_ECAP
    out_varnames.push_back("X(Ne21)");
    out_varnames.push_back("X(F21)");
    out_varnames.push_back("A21_electron_capture_rate");
    out_varnames.push_back("A21_electron_capture_nu_loss");
    nA21=4;
#endif
#endif

//Add A23 nuclei if either is included. 
//Add rates only if react is included
#ifndef SKIP_A23_BETA 
    out_varnames.push_back("X(Na23)");
    out_varnames.push_back("X(Ne23)");
    out_varnames.push_back("A23_beta_decay_rate");   
    out_varnames.push_back("A23_beta_decay_nu_loss");
    nA23=4;
#ifndef SKIP_A23_ECAP
    out_varnames.push_back("A23_electron_capture_rate");
    out_varnames.push_back("A23_electron_capture_nu_loss");
    nA23=6;
#endif
#else
#ifndef SKIP_A23_ECAP
    out_varnames.push_back("X(Na23)");
    out_varnames.push_back("X(Ne23)");
    out_varnames.push_back("A23_electron_capture_rate");
    out_varnames.push_back("A23_electron_capture_nu_loss");
    nA23=4;
#endif
#endif

//Add A23 nuclei if either is included. 
//Add rates only if react is included
#ifndef SKIP_A25_BETA 
    out_varnames.push_back("X(Mg25)");
    out_varnames.push_back("X(Na25)");
    out_varnames.push_back("A25_beta_decay_rate");   
    out_varnames.push_back("A25_beta_decay_nu_loss");
    nA25=4;
#ifndef SKIP_A25_ECAP
    out_varnames.push_back("A25_electron_capture_rate");
    out_varnames.push_back("A25_electron_capture_nu_loss");
    nA25=6;
#endif
#else
#ifndef SKIP_A25_ECAP
    out_varnames.push_back("X(Mg25)");
    out_varnames.push_back("X(Na25)");
    out_varnames.push_back("A25_electron_capture_rate");
    out_varnames.push_back("A25_electron_capture_nu_loss");
    nA25=4;
#endif
#endif

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
     
             
                composition(eos_state);
                eos_to_burn(eos_state, burn_state);
             
                amrex::Real rate, drate_dt, edot_nu, edot_gamma;
                amrex::Real rhoy = burn_state.rho * burn_state.y_e;

                int ine21, if21, ine23, ina23, img25, ina25;

                // do A=21
#ifndef SKIP_A21_ECAP
                if21 = network_spec_index("fluorine-21");
                ine21 = network_spec_index("neon-21");
                tabular_evaluate(j_Ne21_F21_meta, j_Ne21_F21_rhoy, j_Ne21_F21_temp, j_Ne21_F21_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);

                new_arr(i, j, k, 3) = burn_state.xn[ine21] * rate;
                new_arr(i, j, k, 4) = C::Legacy::n_A * burn_state.xn[ine21]/21 * (edot_nu + edot_gamma);
#endif

#ifndef SKIP_A21_BETA 
                if21 = network_spec_index("fluorine-21");
                ine21 = network_spec_index("neon-21");
                tabular_evaluate(j_F21_Ne21_meta, j_F21_Ne21_rhoy, j_F21_Ne21_temp, j_F21_Ne21_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);

                new_arr(i, j, k, 1 + nA21 - 2) = burn_state.xn[if21] * rate;
                new_arr(i, j, k, 1 + nA21 - 1) = C::Legacy::n_A * burn_state.xn[if21]/21 * (edot_nu + edot_gamma);
#endif



                // Do A=23
#ifndef SKIP_A23_ECAP
                ine23 = network_spec_index("neon-23");
                ina23 = network_spec_index("sodium-23");
                tabular_evaluate(j_Na23_Ne23_meta, j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);

                new_arr(i, j, k, nA21+3) = burn_state.xn[ina23] * rate;
                new_arr(i, j, k, nA21+4) = C::Legacy::n_A * burn_state.xn[ina23]/23 * (edot_nu + edot_gamma);
#endif

#ifndef SKIP_A23_BETA 
                ine23 = network_spec_index("neon-23");
                ina23 = network_spec_index("sodium-23");
                tabular_evaluate(j_Ne23_Na23_meta, j_Ne23_Na23_rhoy, j_Ne23_Na23_temp, j_Ne23_Na23_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);
                
                new_arr(i, j, k, 1 + nA21 + nA23 - 2) = burn_state.xn[ine23] * rate;
                new_arr(i, j, k, 1 + nA21 + nA23 - 1) = C::Legacy::n_A * burn_state.xn[ine23]/23 * (edot_nu + edot_gamma);
#endif

                // Do A=25
#ifndef SKIP_A25_ECAP
                ina25 = network_spec_index("sodium-25");
                img25 = network_spec_index("magnesium-25");
                tabular_evaluate(j_Mg25_Na25_meta, j_Mg25_Na25_rhoy, j_Mg25_Na25_temp, j_Mg25_Na25_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);

                new_arr(i, j, k, nA21+nA23+3) = burn_state.xn[img25] * rate;
                new_arr(i, j, k, nA21+nA23+4) = C::Legacy::n_A * burn_state.xn[img25]/25 * (edot_nu + edot_gamma);
#endif

#ifndef SKIP_A25_BETA 
                ina25 = network_spec_index("sodium-25");
                img25 = network_spec_index("magnesium-25");
                tabular_evaluate(j_Na25_Mg25_meta, j_Na25_Mg25_rhoy, j_Na25_Mg25_temp, j_Na25_Mg25_data,
                                  rhoy, burn_state.T, rate, drate_dt, edot_nu, edot_gamma);

                new_arr(i, j, k, 1+nA21+nA23+nA25-2) = burn_state.xn[ina25] * rate;
                new_arr(i, j, k, 1+nA21+nA23+nA25-1) = C::Legacy::n_A * burn_state.xn[ina25]/25 * (edot_nu + edot_gamma);
#endif


                //thermal neutrino loss   
                Real sneut, dsneutdt, dsneutdd, dsnuda, dsnudz;

                constexpr int do_T_derivatives = 0;
                sneut5<do_T_derivatives>(burn_state.T, burn_state.rho, burn_state.abar, burn_state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

                //save values
                new_arr(i, j, k, 0) = rho_arr(i, j, k);
                if (nA21){
                    new_arr(i, j, k, 1) = X_arr(i, j, k, ina23); 
                    new_arr(i, j, k, 2) = X_arr(i, j, k, ine23);
                }
                if (nA23){
                    new_arr(i, j, k, nA21+1) = X_arr(i, j, k, ine21); 
                    new_arr(i, j, k, nA21+2) = X_arr(i, j, k, if21);
                }
                if (nA25){
                    new_arr(i, j, k, nA21+nA23+1) = X_arr(i, j, k, img25); 
                    new_arr(i, j, k, nA21+nA23+2) = X_arr(i, j, k, ina25);
                }


                new_arr(i, j, k, nA21+nA23+nA25) = sneut;
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