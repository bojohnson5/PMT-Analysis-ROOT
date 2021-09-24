#ifndef ROOTER_HPP
#define ROOTER_HPP

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <fstream>

struct Results {
  int num_waveforms;
  double ped_sig;
  double ped_sig_err;
  double pv;
  double res;
  double res_err;
  double adj_spe_peak;
  double adj_spe_peak_err;
  double per_0pe;
};

class Rooter {
private:
  std::string file_name;
  std::string red_f_name;
  std::string plot_f_name;
  TFile *file;
  TTree *tree;

  std::vector<UShort_t> *waveform;
  double baseline;
  double polarity;
  long long timestamp;

  long long tot_ents;
  double time_res;
  int rec_len;
  int run_time;
  bool amp;
  double led;
  bool filter;
  int pmt_v;

  int count_pulses(int thre, int st, int en);
  int count_waveform_pulses(int thre, int st, int en);
  int find_average_peak();

public:
  // constructors and destructors
  Rooter(std::string f_name, int run_time, int rec_len,
         double led, bool filter, int pmt_v, bool amp = false);
  ~Rooter();

  // analysis-related functions
  long long get_tot_ents() { return tot_ents; }
  void view_waveform(int num);
  void view_waveform_raw(int num);
  void view_waveform_wind(int num, double view_wind[2]);
  void view_max_amp(int view_wind[2], int n_bins = 100);
  // void view_spectrum(double wind[2], int view_wind[2], int n_bins = 100);
  void view_spectrum(int view_wind[2], int n_bins = 100);
  void view_multi_spec(std::vector<std::array<int, 2>> winds, int view_wind[2],
                       int n_bins = 100);
  void view_fit_spectrum(double wind[2], int view_wind[2],
                         std::vector<std::array<int, 2>> fits,
                         int n_bins = 100);
  Results fit_spectrum(double wind[2], int view_wind[2],
                       std::vector<std::array<int, 2>> fits, int n_bins = 100);
  void view_dark_rates(int st_thres, int en_thres, int step);
  double dark_rate(int thres);
  void pre_late_pulsing(std::vector<std::array<int, 2>> winds, int spe_thre,
                        int thre, std::ofstream& ofile, std::string name);

  std::string get_red_name() { return red_f_name; }
  bool get_amp() { return amp; }
  void get_timestamps() {
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      std::cout << timestamp << std::endl;
    }
  }
};

#endif // ROOTER_HPP
