/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */


#include <CCA/Components/Examples/ParticleTest1.h>
#include <CCA/Components/Examples/ExamplesLabel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/SimpleMaterial.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/ParticleInterpolator.h>
#include <Core/Grid/LinearInterpolator.h>


#include <iostream>

using namespace std;

using namespace Uintah;

ParticleTest1::ParticleTest1(const ProcessorGroup* myworld)
  : UintahParallelComponent(myworld)
{
  phi_label = VarLabel::create("phi", 
                               NCVariable<double>::getTypeDescription());
  
  residual_label = VarLabel::create("residual", 
                                    sum_vartype::getTypeDescription());
  
  // Create face-centered variables
  FX_label = VarLabel::create("F_x", NCVariable<double>::getTypeDescription());
  
  FY_label = VarLabel::create("F_y", NCVariable<double>::getTypeDescription());
  
  FZ_label = VarLabel::create("F_z", NCVariable<double>::getTypeDescription());
  
  lb_ = scinew ExamplesLabel();
}

ParticleTest1::~ParticleTest1()
{
  VarLabel::destroy(phi_label);
  VarLabel::destroy(residual_label);
  VarLabel::destroy(FX_label);
  VarLabel::destroy(FY_label);
  VarLabel::destroy(FZ_label);
  delete lb_;
}

void ParticleTest1::problemSetup(const ProblemSpecP& params, 
                                 const ProblemSpecP& restart_prob_spec, 
                                 GridP& /*grid*/,SimulationStateP& sharedState)
{
  sharedState_ = sharedState;
  dynamic_cast<Scheduler*>(getPort("scheduler"))->setPositionVar(lb_->pXLabel);
  //dynamic_cast<Scheduler*>(getPort("scheduler"))->setPositionVar(lb_->pVLabel);
  ProblemSpecP pt1 = params->findBlock("ParticleTest1");
  pt1->getWithDefault("doOutput", doOutput_, 0);
  pt1->getWithDefault("doGhostCells", doGhostCells_ , 0);
  ProblemSpecP poisson = params->findBlock("Poisson");
  poisson->require("delt", poisson_delt_);
  poisson->require("maxresidual", poisson_maxresidual_);
  
  mymat_ = scinew SimpleMaterial();
  sharedState_->registerSimpleMaterial(mymat_);

}
 
void ParticleTest1::scheduleInitialize(const LevelP& level,
			       SchedulerP& sched)
{
  Task* task = scinew Task("initialize",
			   this, &ParticleTest1::initialize);
  task->computes(lb_->pXLabel);
  task->computes(lb_->pVxLabel); task->computes(lb_->pVyLabel); task->computes(lb_->pVzLabel);
  task->computes(lb_->pMassLabel);
  task->computes(lb_->pParticleIDLabel);
  task->computes(phi_label);
  task->computes(FX_label);
  task->computes(FY_label);
  task->computes(FZ_label);
  sched->addTask(task, level->eachPatch(), sharedState_->allMaterials());
}
 
void ParticleTest1::scheduleComputeStableTimestep(const LevelP& level,
					  SchedulerP& sched)
{
  Task* task = scinew Task("computeStableTimestep",
			   this, &ParticleTest1::computeStableTimestep);
  task->computes(sharedState_->get_delt_label(),level.get_rep());
  sched->addTask(task, level->eachPatch(), sharedState_->allMaterials());

}

void
ParticleTest1::scheduleTimeAdvance( const LevelP& level, SchedulerP& sched)
{

  const MaterialSet* matls = sharedState_->allMaterials();

  Task* task = scinew Task("timeAdvance",
			   this, &ParticleTest1::timeAdvance, level, sched.get_rep());

  // set this in problemSetup.  0 is no ghost cells, 1 is all with 1 ghost
  // atound-node, and 2 mixes them
  if (doGhostCells_ == 0) {
    task->requires(Task::OldDW, lb_->pParticleIDLabel, Ghost::None, 0);
    task->requires(Task::OldDW, lb_->pXLabel, Ghost::None, 0);
    task->requires(Task::OldDW, lb_->pVxLabel, Ghost::None, 0); task->requires(Task::OldDW, lb_->pVyLabel, Ghost::None, 0); task->requires(Task::OldDW, lb_->pVzLabel, Ghost::None, 0);
    task->requires(Task::OldDW, lb_->pMassLabel, Ghost::None, 0);
  }
  
  else if (doGhostCells_ == 1) {
    task->requires(Task::OldDW, lb_->pXLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, lb_->pVxLabel, Ghost::AroundNodes, 0); task->requires(Task::OldDW, lb_->pVyLabel, Ghost::AroundNodes, 0); task->requires(Task::OldDW, lb_->pVzLabel, Ghost::AroundNodes, 0);
    task->requires(Task::OldDW, lb_->pMassLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, lb_->pParticleIDLabel, Ghost::AroundNodes, 1);
  }
  else if (doGhostCells_ == 2) {
    task->requires(Task::OldDW, lb_->pXLabel, Ghost::None, 0);
    task->requires(Task::OldDW, lb_->pVxLabel, Ghost::None, 0); task->requires(Task::OldDW, lb_->pVyLabel, Ghost::None, 0); task->requires(Task::OldDW, lb_->pVzLabel, Ghost::None, 0);
    task->requires(Task::OldDW, lb_->pMassLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, lb_->pParticleIDLabel, Ghost::None, 0);
  }

  task->computes(lb_->pXLabel_preReloc);
  task->computes(lb_->pVxLabel_preReloc); task->computes(lb_->pVyLabel_preReloc); task->computes(lb_->pVzLabel_preReloc);
  task->computes(lb_->pMassLabel_preReloc);
  task->computes(lb_->pParticleIDLabel_preReloc);
  
  // Gradient calculation (after poisson):
  //task->computes(FX_label);
  //task->computes(FY_label);
  //task->computes(FZ_label);
  
  // Poisson stuff:
  task->hasSubScheduler();
  task->requires(Task::OldDW, phi_label, Ghost::AroundNodes, 1);
  task->computes(phi_label);
  LoadBalancer* lb = sched->getLoadBalancer();
  const PatchSet* perproc_patches = lb->getPerProcessorPatchSet(level);
  sched->addTask(task, perproc_patches, sharedState_->allMaterials());
  
  // Gradient
  //// Calculate gradient:
  Task* gradient_task = scinew Task("calculate_gradient", this, &ParticleTest1::calculate_potential_gradients);
  gradient_task->requires(Task::NewDW, phi_label, Ghost::AroundNodes, 1);
  gradient_task->computes(FX_label);
  gradient_task->computes(FY_label);
  gradient_task->computes(FZ_label);
  sched->addTask(gradient_task, perproc_patches, sharedState_->allMaterials());

  lb_->d_particleState.clear();
  lb_->d_particleState_preReloc.clear();
  for (int m = 0; m < matls->size(); m++) {
    vector<const VarLabel*> vars;
    vector<const VarLabel*> vars_preReloc;

    vars.push_back(lb_->pMassLabel);
    vars.push_back(lb_->pParticleIDLabel);
    vars.push_back(lb_->pVxLabel); vars.push_back(lb_->pVyLabel); vars.push_back(lb_->pVzLabel);

    vars_preReloc.push_back(lb_->pMassLabel_preReloc);
    vars_preReloc.push_back(lb_->pParticleIDLabel_preReloc);
    vars_preReloc.push_back(lb_->pVxLabel_preReloc); vars_preReloc.push_back(lb_->pVyLabel_preReloc); vars_preReloc.push_back(lb_->pVzLabel_preReloc);
    
    lb_->d_particleState.push_back(vars);
    lb_->d_particleState_preReloc.push_back(vars_preReloc);
  }

  sched->scheduleParticleRelocation(level, lb_->pXLabel_preReloc,
				    lb_->d_particleState_preReloc,
				    lb_->pXLabel, lb_->d_particleState,
				    lb_->pParticleIDLabel, matls);

}

void ParticleTest1::computeStableTimestep(const ProcessorGroup* /*pg*/,
				     const PatchSubset* patches,
				     const MaterialSubset* /*matls*/,
				     DataWarehouse*,
				     DataWarehouse* new_dw)
{
  new_dw->put(delt_vartype(1), sharedState_->get_delt_label(),getLevel(patches));
}

void ParticleTest1::initialize(const ProcessorGroup*,
			  const PatchSubset* patches,
			  const MaterialSubset* matls,
			  DataWarehouse* /*old_dw*/, DataWarehouse* new_dw)
{
  for( int p=0; p<patches->size(); ++p ){
    const Patch* patch = patches->get(p);
    const Point low = patch->cellPosition(patch->getCellLowIndex());
    const Point high = patch->cellPosition(patch->getCellHighIndex());
    for(int m = 0;m<matls->size();m++){
      srand(1);
      const int numParticles = 10;
      const int matl = matls->get(m);

      ParticleVariable<Point> px;
      //ParticleVariable<Point> pv;
      ParticleVariable<double> pv[3];
      ParticleVariable<double> pmass;
      ParticleVariable<long64> pids;

      ParticleSubset* subset = new_dw->createParticleSubset(numParticles,matl,patch);
      new_dw->allocateAndPut( px,    lb_->pXLabel,          subset );
      new_dw->allocateAndPut( pv[0], lb_->pVxLabel,         subset );
      new_dw->allocateAndPut( pv[1], lb_->pVyLabel,         subset );
      new_dw->allocateAndPut( pv[2], lb_->pVzLabel,         subset );
      new_dw->allocateAndPut( pmass, lb_->pMassLabel,       subset );
      new_dw->allocateAndPut( pids,  lb_->pParticleIDLabel, subset );

      for( int i = 0; i < numParticles; ++i ){
        const Point pos( (((float) rand()) / RAND_MAX * ( high.x() - low.x()-1) + low.x()),
                         (((float) rand()) / RAND_MAX * ( high.y() - low.y()-1) + low.y()),
                         (((float) rand()) / RAND_MAX * ( high.z() - low.z()-1) + low.z()) );
	//const Point vel( 0.25, 0.0, 0.0 );
	//pv[i] = vel;
        px[i] = pos;
        pv[0][i] = 0.25;
        pv[1][i] = 0;
        pv[2][i] = 0;
        pids[i] = patch->getID()*numParticles+i;
        pmass[i] = ((float) rand()) / RAND_MAX * 10;
      }
    }
  }
  
  // Poisson solver stuff:
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);
      NCVariable<double> phi;
      NCVariable<double> F_x;
      NCVariable<double> F_y;
      NCVariable<double> F_z;
      new_dw->allocateAndPut(F_x, FX_label, matl, patch);
      new_dw->allocateAndPut(F_y, FY_label, matl, patch);
      new_dw->allocateAndPut(F_z, FZ_label, matl, patch);
      new_dw->allocateAndPut(phi, phi_label, matl, patch);
      phi.initialize(0);
      F_x.initialize(0);
      F_y.initialize(0);
      F_z.initialize(0);

      for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
           face=Patch::nextFace(face)) {
        
        if (patch->getBCType(face) == Patch::None) {
          int numChildren = 
            patch->getBCDataArray(face)->getNumberChildren(matl);
          for (int child = 0; child < numChildren; child++) {
            Iterator nbound_ptr, nu;
            
            const BoundCondBase* bcb = patch->getArrayBCValues(face,matl,"Phi",
                                                               nu,nbound_ptr,
                                                               child);
            
            const BoundCond<double>* bc = 
              dynamic_cast<const BoundCond<double>*>(bcb); 
            double value = bc->getValue();
            for (nbound_ptr.reset(); !nbound_ptr.done();nbound_ptr++) {
              phi[*nbound_ptr]=value;
              
            }
          }
        }
      }    
    }
  }
}

void ParticleTest1::calculate_potential_gradients(const ProcessorGroup* pg,
                                                  const PatchSubset* patches,
                                                  const MaterialSubset* matls,
                                                  DataWarehouse* old_dw, DataWarehouse* new_dw) {

  //cout << "hello" << endl;
  // Poisson solver stuff:
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);

      constNCVariable<double> phi;
      NCVariable<double> F_x, F_y, F_z;
      
      new_dw->get(phi, phi_label, matl, patch, Ghost::AroundNodes, 1);
      new_dw->allocateAndPut(F_x, FX_label, matl, patch);
      new_dw->allocateAndPut(F_y, FY_label, matl, patch);
      new_dw->allocateAndPut(F_z, FZ_label, matl, patch);
      F_x.initialize(0);
      F_y.initialize(0);
      F_z.initialize(0);
      
      // Calculate gradient:
      Vector dx = patch->dCell();
      IntVector l = patch->getNodeLowIndex();
      IntVector h = patch->getNodeHighIndex();
      
      l += IntVector(patch->getBCType(Patch::xminus) == Patch::Neighbor?0:1,
                     patch->getBCType(Patch::yminus) == Patch::Neighbor?0:1,
                     patch->getBCType(Patch::zminus) == Patch::Neighbor?0:1);
      h -= IntVector(patch->getBCType(Patch::xplus) == Patch::Neighbor?0:1,
                     patch->getBCType(Patch::yplus) == Patch::Neighbor?0:1,
                     patch->getBCType(Patch::zplus) == Patch::Neighbor?0:1);

      for( NodeIterator iter(l,h); !iter.done(); iter++ ) {
        F_x[*iter] = (phi[*iter+IntVector(1,0,0)] - phi[*iter + IntVector(-1,0,0)]) / (2.0*dx[0]);
        F_y[*iter] = (phi[*iter+IntVector(0,1,0)] - phi[*iter + IntVector(0,-1,0)]) / (2.0*dx[1]);
        F_z[*iter] = (phi[*iter+IntVector(0,0,1)] - phi[*iter + IntVector(0,0,-1)]) / (2.0*dx[2]);
      }
    }
  }
//  for(int p=0;p<patches->size();p++){
//    const Patch* patch = patches->get(p);
//    for(int m = 0;m<matls->size();m++){
//      int matl = matls->get(m);
//      constNCVariable<double> phi;
//      old_dw->get(phi, phi_label, matl, patch, Ghost::AroundNodes, 1);
//      NCVariable<double> newphi;
//      new_dw->allocateAndPut(newphi, phi_label, matl, patch);
//      newphi.copyPatch(phi, newphi.getLow(), newphi.getHigh());
//      double residual=0;
//      IntVector l = patch->getNodeLowIndex();
//      IntVector h = patch->getNodeHighIndex();
//      l += IntVector(patch->getBCType(Patch::xminus) == Patch::Neighbor?0:1,
//		     patch->getBCType(Patch::yminus) == Patch::Neighbor?0:1,
//		     patch->getBCType(Patch::zminus) == Patch::Neighbor?0:1);
//      h -= IntVector(patch->getBCType(Patch::xplus) == Patch::Neighbor?0:1,
//		     patch->getBCType(Patch::yplus) == Patch::Neighbor?0:1,
//		     patch->getBCType(Patch::zplus) == Patch::Neighbor?0:1);
//      for(NodeIterator iter(l, h);!iter.done(); iter++){
//	newphi[*iter]=(1./6)*(
//	  phi[*iter+IntVector(1,0,0)]+phi[*iter+IntVector(-1,0,0)]+
//	  phi[*iter+IntVector(0,1,0)]+phi[*iter+IntVector(0,-1,0)]+
//	  phi[*iter+IntVector(0,0,1)]+phi[*iter+IntVector(0,0,-1)]);
//	double diff = newphi[*iter]-phi[*iter];
//	residual += diff*diff;
//      }
//      new_dw->put(sum_vartype(residual), residual_label);
//    }
//  }
//template <class T>
//void compute_Mag_gradient( constCCVariable<T>& q_CC,
//                           CCVariable<T>& mag_grad_q_CC,
//                           const Patch* patch)
//{
//  Vector dx = patch->dCell();
//  for(CellIterator iter = patch->getCellIterator();!iter.done();iter++){
//    IntVector c = *iter;
//    Vector grad_q_CC;
//
//    for(int dir = 0; dir <3; dir ++ ) {
//      IntVector r = c;
//      IntVector l = c;
//      T inv_dx = (T) 0.5 /dx[dir];
//      r[dir] += 1;
//      l[dir] -= 1;
//      grad_q_CC[dir] = (q_CC[r] - q_CC[l])*inv_dx;
//    }
//    mag_grad_q_CC[c] = grad_q_CC.length();
//  }
  //old_dw->get
  //Interpolate Particles To Grid (from AMRMPM)
//      //double pSp_vol = 1./mpm_matl->getInitialDensity();
//      for (ParticleSubset::iterator iter = pset->begin();iter != pset->end(); iter++){
//        particleIndex idx = *iter;
//
//        // Get the node indices that surround the cell
//        interpolator->findCellAndWeights(px[idx],ni,S,psize[idx],pDeformationMeasure[idx]); //ParticleInterpolator
//
//        pmom = pvelocity[idx]*pmass[idx];
//        
//        // Add each particles contribution to the local mass & velocity 
//        IntVector node;
//        for(int k = 0; k < n8or27; k++) {
//          node = ni[k];
//          if(patch->containsNode(node)) {
//            gmass[node]          += pmass[idx]                     * S[k];
//            gvelocity[node]      += pmom                           * S[k];
//            gvolume[node]        += pvolume[idx]                   * S[k];
//            gexternalforce[node] += pexternalforce[idx]            * S[k];
//            gTemperature[node]   += pTemperature[idx] * pmass[idx] * S[k];
//          }
//        }
//      }  // End of particle loop
}


void ParticleTest1::timeAdvance(const ProcessorGroup* pg,
			   const PatchSubset* patches,
			   const MaterialSubset* matls,
			   DataWarehouse* old_dw, DataWarehouse* new_dw,
			   LevelP level, Scheduler* sched)
{
  SchedulerP subsched = sched->createSubScheduler();
  subsched->initialize();
  GridP grid = level->getGrid();

  // Create the tasks
  Task* task = scinew Task("solve_poisson",
			   this, &ParticleTest1::poisson_solver);
  task->requires(Task::OldDW, phi_label, Ghost::AroundNodes, 1);
  task->computes(phi_label);
  task->computes(residual_label);
  

  
  subsched->addTask(task, level->eachPatch(), sharedState_->allMaterials());
  
  // Compile the scheduler
  subsched->advanceDataWarehouse(grid);
  subsched->compile();
  
  int count = 0;
  double residual;
  subsched->get_dw(1)->transferFrom(old_dw, phi_label, patches, matls);
  // Loop poisson iterations
  do {
    subsched->advanceDataWarehouse(grid);
    subsched->get_dw(0)->setScrubbing(DataWarehouse::ScrubComplete);
    subsched->get_dw(1)->setScrubbing(DataWarehouse::ScrubNonPermanent);
    subsched->execute();    
  
    sum_vartype residual_var;
    subsched->get_dw(1)->get(residual_var, residual_label);
    residual = residual_var;
  
    if(pg->myrank() == 0)
      cerr << "Iteration " << count++ << ", residual=" << residual << '\n';
  } while(residual > poisson_maxresidual_);
  
  new_dw->transferFrom(subsched->get_dw(1), phi_label, patches, matls);
  //old_dw->transferFrom(subsched->get_dw(1), phi_label, patches, matls);
  
  //SchedulerP gradient_subsched = sched->createSubScheduler();
  //gradient_subsched->initialize();
  //gradient_subsched->addTask(gradient_task, level->eachPatch(), sharedState_->allMaterials());
  //// Compile the scheduler
  //gradient_subsched->advanceDataWarehouse(grid);
  //gradient_subsched->compile();
  //gradient_subsched->get_dw(1)->transferFrom(old_dw, phi_label, patches, matls);
  //gradient_subsched->advanceDataWarehouse(grid);
  //gradient_subsched->get_dw(0)->setScrubbing(DataWarehouse::ScrubComplete);
  //gradient_subsched->get_dw(1)->setScrubbing(DataWarehouse::ScrubNonPermanent);
  //gradient_subsched->execute();
  //new_dw->transferFrom(gradient_subsched->get_dw(1), FX_label, patches, matls);
  //new_dw->transferFrom(gradient_subsched->get_dw(1), FY_label, patches, matls);
  //new_dw->transferFrom(gradient_subsched->get_dw(1), FZ_label, patches, matls);
  
  // Now Poisson has been calculated, calculate the gradients
  //calculate_potential_gradients(pg, patches, matls, old_dw, new_dw);
  
  // Advance particles

  for( int p=0; p<patches->size(); ++p ){
    const Patch* patch = patches->get(p);
    ParticleInterpolator* interpolator=scinew LinearInterpolator(patch);
    for( int m = 0; m<matls->size(); ++m ){
      int matl = matls->get(m);
      ParticleSubset* pset = old_dw->getParticleSubset(matl, patch);
      ParticleSubset* delset = scinew ParticleSubset(0,matl,patch);

      // Get the arrays of particle values to be changed
      constParticleVariable<Point> px;
      ParticleVariable<Point> pxnew;
      //constParticleVariable<Point> pv;
      //ParticleVariable<Point> pvnew;
      constParticleVariable<double> pv[3];
      ParticleVariable<double> pvnew[3];
      constParticleVariable<double> pmass;
      ParticleVariable<double> pmassnew;
      constParticleVariable<long64> pids;
      ParticleVariable<long64> pidsnew;

      old_dw->get(pmass, lb_->pMassLabel,               pset);
      old_dw->get(px,    lb_->pXLabel,                  pset);
      old_dw->get(pv[0], lb_->pVxLabel,                 pset);
      old_dw->get(pv[1], lb_->pVyLabel,                 pset);
      old_dw->get(pv[2], lb_->pVzLabel,                 pset);
      old_dw->get(pids,  lb_->pParticleIDLabel,         pset);

      new_dw->allocateAndPut(pmassnew, lb_->pMassLabel_preReloc,       pset);
      new_dw->allocateAndPut(pxnew,    lb_->pXLabel_preReloc,          pset);
      new_dw->allocateAndPut(pvnew[0], lb_->pVxLabel_preReloc,         pset);
      new_dw->allocateAndPut(pvnew[1], lb_->pVyLabel_preReloc,         pset);
      new_dw->allocateAndPut(pvnew[2], lb_->pVzLabel_preReloc,         pset);
      
      new_dw->allocateAndPut(pidsnew,  lb_->pParticleIDLabel_preReloc, pset);

      // Update particle position and velocity:
      for( unsigned i = 0; i < pset->numParticles(); ++i ){
//        particleIndex idx = i;
//        // Calculate gradient of potential:
//        {
//          Vector dx = patch->dCell();
//          double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
//          vector<IntVector> ni(interpolator->size());
//          vector<double> S(interpolator->size());
//          vector<Vector> d_S(interpolator->size());
//          // Define gradient
//          Matrix3 gradient_of_potential(0.0);
//          // Get nodes
//          vector<IntVector> node_indices; 
//          interpolator->findCellAndShapeDerivatives(px[idx], node_indices, d_S, psize[idx], pFOld[idx]);
//        }
        
        //Point pos( px[i].x() + pv[i].x(), px[i].y() + pv[i].y(), px[i].z() + pv[i].z());
        Point pos( px[i].x() + pv[0][i], px[i].y() + pv[1][i], px[i].z() + pv[2][i]);
        pxnew[i] = pos;
        for( int component = 0; component < 3; ++component ) {
          pvnew[component][i] = pv[component][i];
        }
        pidsnew[i] = pids[i];
        pmassnew[i] = pmass[i];
        if (doOutput_)
          cout << " Patch " << patch->getID() << ": ID " 
               << pidsnew[i] << ", pos " << pxnew[i] 
               << ", mass " << pmassnew[i] << endl;
      }
      new_dw->deleteParticles(delset);
    }
    delete interpolator;
  }
}


void ParticleTest1::poisson_solver(const ProcessorGroup*,
		       const PatchSubset* patches,
		       const MaterialSubset* matls,
		       DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);
      constNCVariable<double> phi;
      old_dw->get(phi, phi_label, matl, patch, Ghost::AroundNodes, 1);
      NCVariable<double> newphi;
      new_dw->allocateAndPut(newphi, phi_label, matl, patch);
      newphi.copyPatch(phi, newphi.getLow(), newphi.getHigh());
      double residual=0;
      IntVector l = patch->getNodeLowIndex();
      IntVector h = patch->getNodeHighIndex();
      l += IntVector(patch->getBCType(Patch::xminus) == Patch::Neighbor?0:1,
		     patch->getBCType(Patch::yminus) == Patch::Neighbor?0:1,
		     patch->getBCType(Patch::zminus) == Patch::Neighbor?0:1);
      h -= IntVector(patch->getBCType(Patch::xplus) == Patch::Neighbor?0:1,
		     patch->getBCType(Patch::yplus) == Patch::Neighbor?0:1,
		     patch->getBCType(Patch::zplus) == Patch::Neighbor?0:1);
      for(NodeIterator iter(l, h);!iter.done(); iter++){
	newphi[*iter]=(1./6)*(
	  phi[*iter+IntVector(1,0,0)]+phi[*iter+IntVector(-1,0,0)]+
	  phi[*iter+IntVector(0,1,0)]+phi[*iter+IntVector(0,-1,0)]+
	  phi[*iter+IntVector(0,0,1)]+phi[*iter+IntVector(0,0,-1)]);
	double diff = newphi[*iter]-phi[*iter];
	residual += diff*diff;
      }
      new_dw->put(sum_vartype(residual), residual_label);
    }
  }
}