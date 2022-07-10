//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>

#include <AnalysisTree/AnalysisTask.hpp>

#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Particle.hpp>
#include <AnalysisTree/BranchConfig.hpp>
#include <AnalysisTree/EventHeader.hpp>

#include <memory>
#include <string>
#include <TProfile.h>
#include <TProfile2D.h>

namespace AnalysisTree {

struct qvector{
  double x;
  double y;
  double m;

  double operator*( qvector other) const{
    return (x*other.x + y*other.y) / (m*other.m);
  }
};

class TracksProcessor : public Task {

public:
  void Init() override;
  void Exec() override;
  void Finish() override;
  void SetOutFileName(const std::string &out_file_name) {
    out_file_name_ = out_file_name;
  }

protected:
  void LoopSimParticles();
  void FillReferenceVectors();

  Branch in_particles_;
  Branch in_event_header_;

  std::vector<double> b_edges_;

  double centrality_{-1};

  qvector psi_vector_;
  qvector f1_vector_;
  qvector f2_vector_;
  qvector f3_vector_;

  TProfile* p1_f1_psi_;
  TProfile* p1_f2_psi_;
  TProfile* p1_f3_psi_;

  TProfile* p1_f1_f2_;
  TProfile* p1_f1_f3_;
  TProfile* p1_f2_f3_;

  TProfile2D*p2_proton_psi_;
  TProfile2D*p2_proton_f1_;
  TProfile2D*p2_proton_f2_;
  TProfile2D*p2_proton_f3_;

  TH1F* h1_eta_;
  TH1F* h1_centrality_;

  std::string out_file_name_;
  TFile* out_file_;

};

}
#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
