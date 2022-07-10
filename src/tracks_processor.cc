//
// Created by mikhail on 11/30/20.
//

#include "tracks_processor.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TaskManager.hpp"

#include <AnalysisTree/DataHeader.hpp>
#include <random>
#include <TH2F.h>

namespace AnalysisTree {

void TracksProcessor::Init() {
  auto man = TaskManager::GetInstance();
  auto chain = man->GetChain();
  auto data_header = man->GetDataHeader();

  AddInputBranch("particles");
  AddInputBranch("event_header");

  in_particles_ = chain->GetBranch("particles");
  in_particles_.Freeze();
  in_event_header_ = chain->GetBranch("event_header");
  in_event_header_.Freeze();

  p1_f1_psi_ = new TProfile( "p1_f1_psi_", ";centrality (%)", 10, 0, 100 );
  p1_f2_psi_ = new TProfile( "p1_f2_psi_", ";centrality (%)", 10, 0, 100 );
  p1_f3_psi_ = new TProfile( "p1_f3_psi_", ";centrality (%)", 10, 0, 100 );

  p1_f1_f2_ = new TProfile( "p1_f1_f2_", ";centrality (%)", 10, 0, 100 );
  p1_f1_f3_ = new TProfile( "p1_f1_f3_", ";centrality (%)", 10, 0, 100 );
  p1_f2_f3_ = new TProfile( "p1_f2_f3_", ";centrality (%)", 10, 0, 100 );

  p2_proton_psi_ = new TProfile2D( "p2_proton_psi_", ";centrality (%);y_{cm}",
                                 10, 0, 100,
                                 20, -1, 1);
  p2_proton_f1_ = new TProfile2D( "p2_proton_f1_", ";centrality (%);y_{cm}",
                                 10, 0, 100,
                                 20, -1, 1);
  p2_proton_f2_ = new TProfile2D( "p2_proton_f2_", ";centrality (%);y_{cm}",
                                 10, 0, 100,
                                 20, -1, 1);
  p2_proton_f3_ = new TProfile2D( "p2_proton_f3_", ";centrality (%);y_{cm}",
                                 10, 0, 100,
                                 20, -1, 1);

  h1_eta_ = new TH1F( "h1_eta_", ";eta", 350, -1, 6 );
  h1_centrality_ = new TH1F( "h1_centrality_", ";centrality", 10, 0, 100 );

  auto colliding_system = data_header->GetSystem();
  if(colliding_system == "Au+Au" ){
    b_edges_ = {0, 3.888, 5.67, 6.966, 8.1, 9.072, 10.044, 10.854, 11.664, 12.474, 16.2 };
  } else if(colliding_system == "Xe+Cs" ){
    b_edges_ = { 0, 3.608, 5.248, 6.56, 7.708, 8.692, 9.512, 10.332, 11.152, 12.3, 16.4 };
  } else if(colliding_system == "Ag+Ag" ){
    b_edges_ = {  0, 3.12, 4.55, 5.59, 6.5, 7.28, 8.06, 8.71, 9.36, 10.14, 13  };
  }
}

void TracksProcessor::Exec() {
  auto field_centrality = in_event_header_.GetField("centrality");
  auto field_psi_rp = in_event_header_.GetField("psi_rp");
  auto psi_rp = in_event_header_.GetDataRaw<EventHeader*>()->GetField<float>(field_psi_rp.GetFieldId());

  centrality_ = in_event_header_.GetDataRaw<EventHeader*>()->GetField<float>(field_centrality.GetFieldId());
  psi_vector_ = { cos(psi_rp), sin(psi_rp), 1 };

  this->FillReferenceVectors();

  h1_centrality_->Fill(centrality_);
  if( fabs(f1_vector_.m) > std::numeric_limits<float>::min() )
    p1_f1_psi_->Fill( centrality_, f1_vector_*psi_vector_ );
  if( fabs(f2_vector_.m) > std::numeric_limits<float>::min() )
    p1_f2_psi_->Fill( centrality_, f2_vector_*psi_vector_ );
  if( fabs(f3_vector_.m) > std::numeric_limits<float>::min() )
    p1_f3_psi_->Fill( centrality_, f3_vector_*psi_vector_ );

  if( fabs(f1_vector_.m) > std::numeric_limits<float>::min() &&
      fabs(f2_vector_.m) > std::numeric_limits<float>::min())
    p1_f1_f2_->Fill( centrality_, f1_vector_*f2_vector_ );
  if( fabs(f1_vector_.m) > std::numeric_limits<float>::min() &&
      fabs(f3_vector_.m) > std::numeric_limits<float>::min())
    p1_f1_f3_->Fill( centrality_, f1_vector_*f3_vector_ );
  if( fabs(f2_vector_.m) > std::numeric_limits<float>::min() &&
      fabs(f3_vector_.m) > std::numeric_limits<float>::min())
    p1_f2_f3_->Fill( centrality_, f2_vector_*f3_vector_ );

  this->LoopSimParticles();
}

void TracksProcessor::LoopSimParticles() {
  auto field_ycm = in_particles_.GetField("y_cm");
  auto field_phi = in_particles_.GetField("phi");
  auto field_pid = in_particles_.GetField("pid");

  for (size_t i=0; i< in_particles_.size(); ++i) {
    auto in_particle = in_particles_[i];
    auto y_cm = in_particle[field_ycm];
    auto phi = in_particle[field_phi];
    auto pid = int(in_particle[field_pid]);

    if( pid != 2212 )
      continue;

    qvector u_vector{ cos(phi), sin(phi), 1 };

    p2_proton_psi_->Fill( centrality_, y_cm, u_vector*psi_vector_ );
    if( fabs(f1_vector_.m) > std::numeric_limits<float>::min() )
      p2_proton_f1_->Fill( centrality_, y_cm, u_vector*f1_vector_ );
    if( fabs(f2_vector_.m) > std::numeric_limits<float>::min() )
      p2_proton_f2_->Fill( centrality_, y_cm, u_vector*f2_vector_ );
    if( fabs(f3_vector_.m) > std::numeric_limits<float>::min() )
      p2_proton_f3_->Fill( centrality_, y_cm, u_vector*f3_vector_ );
  }
}

void TracksProcessor::FillReferenceVectors() {
  auto field_eta = in_particles_.GetField("eta");
  auto field_phi = in_particles_.GetField("phi");
  auto field_Ekin = in_particles_.GetField("Ekin");

  f1_vector_ = {0.0, 0.0, 0.0};
  f2_vector_ = {0.0, 0.0, 0.0};
  f3_vector_ = {0.0, 0.0, 0.0};

  for (size_t i=0; i< in_particles_.size(); ++i) {
    auto in_particle = in_particles_[i];
    auto eta = in_particle[field_eta];
    auto phi = in_particle[field_phi];
    auto Ekin = in_particle[field_Ekin];

    h1_eta_->Fill(eta);

    if( 4.4 < eta && eta < 5.5 ) {
      f1_vector_.x += Ekin*cos(phi);
      f1_vector_.y += Ekin*sin(phi);
      f1_vector_.m += Ekin;
    }
    if( 3.9 < eta && eta < 4.4 ) {
      f2_vector_.x += Ekin*cos(phi);
      f2_vector_.y += Ekin*sin(phi);
      f2_vector_.m += Ekin;
    }
    if( 3.1 < eta && eta < 3.9 ) {
      f3_vector_.x += Ekin*cos(phi);
      f3_vector_.y += Ekin*sin(phi);
      f3_vector_.m += Ekin;
    }
  }

}
void TracksProcessor::Finish() {
  out_file_ = TFile::Open( out_file_name_.c_str(), "recreate" );

  h1_centrality_->Write();
  h1_eta_->Write();

  p1_f1_psi_->Write();
  p1_f2_psi_->Write();
  p1_f3_psi_->Write();

  p1_f1_f2_->Write();
  p1_f1_f3_->Write();
  p1_f2_f3_->Write();

  p2_proton_psi_->Write();
  p2_proton_f1_->Write();
  p2_proton_f2_->Write();
  p2_proton_f3_->Write();

  out_file_->Close();
}
} // namespace AnalysisTree