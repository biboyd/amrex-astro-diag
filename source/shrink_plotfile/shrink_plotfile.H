//
// This program takes a 3-d cartesian plotfile and calculates
// reduces teh number of saved variables
//
// 
#ifndef shrink_pfile_H
#define shrink_pfile_H

#include <sstream>

#include <extern_parameters.H>
#include <eos.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <string>

void shrink_pfile(){

   //get plotdata   
   PlotFileData pf(plotfilename);
   const int nlevs = pf.finestLevel() + 1;

   // init vectors needed to writeout
   Vector<std::string> var_names = pf.varNames();
   Vector<int> lev_steps(nlevs);
   Vector<IntVect> ref_ratio(nlevs);
   Vector<MultiFab> new_mfs(nlevs);
   Vector<Geometry> geoms(nlevs);

   // add wanted names
   Vector<std::string> small_var_names;
   Vector<int> small_var_idx;
   small_var_names.push_back("rho");
   small_var_names.push_back("entropy");
   small_var_names.push_back("tfromp");
   small_var_names.push_back("p0pluspi");
   small_var_names.push_back("vort");
   small_var_names.push_back("velx");
   small_var_names.push_back("vely");
   small_var_names.push_back("velz");
   small_var_names.push_back("Hnuc");
   small_var_names.push_back("X(c12)");
   small_var_names.push_back("X(h1)");
   small_var_names.push_back("X(he4)");
   small_var_names.push_back("X(mg23)");
   small_var_names.push_back("X(n)");
   small_var_names.push_back("X(na23)");
   small_var_names.push_back("X(ne20)");
   small_var_names.push_back("X(ne23)");
   small_var_names.push_back("X(o16)");
   small_var_names.push_back("omegadot(na23)");
   small_var_names.push_back("omegadot(ne23)");
   small_var_names.push_back("omegadot(c12)");
  
   // loop over smll varnames
   for (auto small_it = small_var_names.begin(); small_it < small_var_names.end(); ++small_it ){
      //loop over all variables to get idx
      auto idx = std::find(var_names.cbegin(), var_names.cend(), *small_it);
      if (idx == var_names.cend()) {
          amrex::Print() << *small_it << std::endl;
          amrex::Error("Error: could not find a variable");
      }
      int var_idx = std::distance(var_names.cbegin(), idx);
      small_var_idx.push_back( var_idx );

   }

   for (int ilev = pf.finestLevel(); ilev >= 0; --ilev) {
      // read plotfile data
      lev_steps[ilev] = pf.levelStep(ilev);
      ref_ratio[ilev] = IntVect(pf.refRatio(ilev), 
                                pf.refRatio(ilev), 
                                pf.refRatio(ilev));

      const MultiFab mf = pf.get(ilev);

      // define new mf with as many componenets as vector
      new_mfs[ilev].define(mf.boxArray(),
                      mf.DistributionMap(),
                      small_var_names.size(),
                      mf.nGrow());

      // copy data into new_mfs[ilev]
      //MultiFab::Copy(new_mfs[ilev], mf, 0, 0, mf.nComp(), mf.nGrow());
      
      // construct geometry from plotfile
      const amrex::RealBox rb(pf.probLo(), pf.probHi());
      const Array<int,3> is_per({0, 0, 0});
      geoms[ilev].define(pf.probDomain(ilev), rb, pf.coordSys(), is_per);

      for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi){
         // find bounds
         const auto tileBox = mfi.tilebox();
         const auto lo = amrex::lbound(tileBox);
         const auto hi = amrex::ubound(tileBox);
         
         // load arrays
         Array4<const Real> const& mf_arr = mf.array(mfi);
         Array4<Real> const& small_arr = new_mfs[ilev].array(mfi);

         // loop over all cells
         ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

            for (auto n=0; n < small_var_names.size(); ++n){
               small_arr(i, j, k, n) = mf_arr(i, j, k, small_var_idx[n]);
            }
         });
      }            
   }
   //writeout data stuf
   const Vector<Geometry> write_geoms = geoms;

   //write out file
   if (writeplot){
      // check if single level
      if (nlevs == 1){

         WriteSingleLevelPlotfile(outfile, 
                                 new_mfs[0], 
                                 small_var_names, 
                                 geoms[0], 
                                 pf.time(),
                                 lev_steps[0]);
      }
      else{
         // write out with gradient fields
         WriteMultiLevelPlotfile(outfile,
                                 nlevs,
                                 GetVecOfConstPtrs(new_mfs),
                                 small_var_names,
                                 write_geoms,
                                 pf.time(),
                                 lev_steps,
                                 ref_ratio);

      }
   }
}
#endif
