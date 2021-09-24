#include "rooter.hpp"
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <sstream>

// bool amp is used if an amplifier was used on the run or not
Rooter::Rooter(std::string f_name, int run_time, int rec_len,
               double led, bool filter, int pmt_v, bool amp)
    : waveform(nullptr), baseline(0.0), polarity(0.0),
      time_res(4.0), run_time(run_time), rec_len(rec_len),
      file_name(f_name), led(led), filter(filter),
      pmt_v(pmt_v), amp(amp) {
  plot_f_name = file_name.substr(file_name.length() - 7, 2);
  red_f_name = file_name.substr(file_name.length() - 7, 7);
  file = new TFile(f_name.c_str());
  tree = static_cast<TTree *>(file->Get("waveformTree"));
  tree->SetBranchAddress("waveform", &waveform);
  tree->SetBranchAddress("baseline", &baseline);
  tree->SetBranchAddress("polarity", &polarity);
  tree->SetBranchAddress("timestamp", &timestamp);
  tot_ents = tree->GetEntries();
}

Rooter::~Rooter() { delete file; }

double Rooter::dark_rate(int thres) {
  int pulses = count_pulses(thres, 0, rec_len);
  double rate = static_cast<double>(pulses) / (run_time * 60);
  return rate;
}

void Rooter::view_dark_rates(int st_thres, int en_thres, int step) {
  int num_thres = (en_thres - st_thres) / step;
  double rates[num_thres];
  double thres[num_thres];
  for (int thre = st_thres, i = 0; thre < en_thres; thre += step, i++) {
    int pulses = count_pulses(thre, 0, rec_len);
    rates[i] = static_cast<double>(pulses) / (run_time * 60);
    thres[i] = thre;
  }
  TCanvas *c = new TCanvas("dark", "Dark Rates", 800, 600);
  TGraph *gr = new TGraph(num_thres, thres, rates);
  gr->GetXaxis()->SetTitle("Thresholds [ADC]");
  gr->GetYaxis()->SetTitle("Rate [Hz]");
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Dark Rates";
  gr->SetTitle(ss.str().c_str());
  gr->Draw("AL*");
  c->Update();
  c->Draw();
}

int Rooter::count_pulses(int thre, int st, int en) {
  int count = 0;
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (st < 0)
      st = 0;
    if (en > waveform->size())
      en = waveform->size();
    bool pulse = false;
    for (int j = st; j < en; j++) {
      int adc = (waveform->at(j) - baseline) * polarity;
      if (!pulse && adc > thre) {
        pulse = true;
        count++;
      }
      if (pulse && adc < thre) {
        pulse = false;
      }
    }
  }
  return count;
}

int Rooter::count_waveform_pulses(int thre, int st, int en) {
  int count = 0;
  if (st < 0)
    st = 0;
  if (en > waveform->size())
    en = waveform->size();
  bool pulse = false;
  for (int i = st; i < en; i++) {
    int adc = (waveform->at(i) - baseline) * polarity;
    if (!pulse && adc > thre) {
      pulse = true;
      count++;
    }
    if (pulse && adc < thre) {
      pulse = false;
    }
  }
  return count;
}

// void Rooter::view_spectrum(double wind[2], int view_wind[2], int n_bins) {
//   int st = wind[0] / time_res;
//   int en = wind[1] / time_res;
//   TCanvas *c = new TCanvas("spec", "Spectrum", 800, 600);
//   c->SetLogy();
//   std::stringstream ss;
//   ss << "Run " << plot_f_name << " Spectrum";
//   TH1D *h = new TH1D("spec_h", ss.str().c_str(), n_bins, view_wind[0], view_wind[1]);
//   h->GetYaxis()->SetTitle("Counts");
//   h->GetXaxis()->SetTitle("Integrated ADC");
//   int total_entries = tree->GetEntries();
//   double integral;
//   for (int i = 0; i < total_entries; i++) {
//     integral = 0.0;
//     tree->GetEntry(i);
//     for (int j = st; j < en + 1; j++) {
//       integral += (waveform->at(j) - baseline) * polarity;
//     }
//     h->Fill(integral);
//   }
//   h->Draw();
//   c->Draw();
//   c->Update();
//   gSystem->ProcessEvents();
// }

void Rooter::view_spectrum(int view_wind[2], int n_bins) {
  int avg_peak = find_average_peak();
  int st = (avg_peak - 12) / time_res;
  int en = (avg_peak + 18) / time_res;
  TCanvas *c = new TCanvas("spec", "Spectrum", 800, 600);
  c->SetLogy();
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Spectrum";
  TH1D *h = new TH1D("spec_h", ss.str().c_str(), n_bins, view_wind[0], view_wind[1]);
  h->GetYaxis()->SetTitle("Counts");
  h->GetXaxis()->SetTitle("Integrated ADC");
  int total_entries = tree->GetEntries();
  double integral;
  for (int i = 0; i < total_entries; i++) {
    integral = 0.0;
    tree->GetEntry(i);
    for (int j = st; j < en + 1; j++) {
      integral += (waveform->at(j) - baseline) * polarity;
    }
    h->Fill(integral);
  }
  h->Draw();
  c->Draw();
  c->Update();
  gSystem->ProcessEvents();
}

void Rooter::view_multi_spec(std::vector<std::array<int, 2>> winds,
                             int view_wind[2], int n_bins) {
  TCanvas *c = new TCanvas("multispec", "Multiple Spectra", 800, 600);
  c->SetLogy();
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Spectra";
  int st, en;
  int total_entries = tree->GetEntries();
  TLegend *l = new TLegend(0.1, 0.8, 0.2, 0.9);
  std::stringstream label;
  for (int k = 0; k < winds.size(); k++) {
    label << "multispec_h" << k;
    TH1D *h =
      new TH1D(label.str().c_str(), ss.str().c_str(), n_bins, view_wind[0], view_wind[1]);
    h->SetStats(kFALSE);
    label.str("");
    ss.str("");
    ss << winds[k][0] << "-" << winds[k][1];
    l->AddEntry(h, ss.str().c_str(), "l");
    h->SetLineColor(k + 1);
    h->GetYaxis()->SetTitle("Counts");
    h->GetXaxis()->SetTitle("Integrated ADC");
    st = winds[k][0] / time_res;
    en = winds[k][1] / time_res;
    double integral;
    for (int i = 0; i < total_entries; i++) {
      integral = 0.0;
      tree->GetEntry(i);
      for (int j = st; j < en + 1; j++) {
        integral += (waveform->at(j) - baseline) * polarity;
      }
      h->Fill(integral);
    }
    if (k == 0) {
      h->Draw();
    } else {
      h->Draw("same");
    }
  }
  l->Draw();
  c->Draw();
  c->Update();
  gSystem->ProcessEvents();
}

void Rooter::view_fit_spectrum(double wind[2], int view_wind[2],
                               std::vector<std::array<int, 2>> fits,
                               int n_bins) {
  int st = wind[0] / time_res;
  int en = wind[1] / time_res;
  TCanvas *c = new TCanvas("spec", "Spectrum", 800, 600);
  c->SetLogy();
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Spectrum";
  TH1D *h = new TH1D("spec_h", ss.str().c_str(), n_bins, view_wind[0], view_wind[1]);
  h->GetYaxis()->SetTitle("Counts");
  h->GetXaxis()->SetTitle("Integrated ADC");
  int total_entries = tree->GetEntries();
  double integral;
  for (int i = 0; i < total_entries; i++) {
    integral = 0.0;
    tree->GetEntry(i);
    for (int j = st; j < en + 1; j++) {
      integral += (waveform->at(j) - baseline) * polarity;
    }
    h->Fill(integral);
  }
  if (fits.size() == 2) {
    TF1 *fit0 = new TF1("fit0", "gaus", fits[0][0], fits[0][1]);
    TF1 *fit1 = new TF1("fit1", "gaus", fits[1][0], fits[1][1]);
    /* guas function of the form f(x)=p0*exp(-0.5*((x-p1)/p2)^2) */
    h->Fit(fit0, "RQ0");
    h->Fit(fit1, "RQ0+");
    double p00 = fit0->GetParameter(0);
    double p01 = fit0->GetParameter(1);
    double p02 = fit0->GetParameter(2);
    double p10 = fit1->GetParameter(0);
    double p11 = fit1->GetParameter(1);
    double p12 = fit1->GetParameter(2);
    double adj_spe_peak = p11 - p01;
    double res = p12 / adj_spe_peak;
    /* find peak-to-valley */
    h->GetXaxis()->SetRangeUser(p01, p11);
    double min = h->GetMinimum();
    int min_bin = h->GetMinimumBin();
    h->GetXaxis()->SetRange(h->GetMinimumBin(), 5000);
    double max = h->GetMaximum();
    double pv = max / min;
    /* find peak-to-valley */
    /* find % of 0-PEs */
    double sum = 0.0;
    for (int i = 0; i < min_bin + 1; i++) {
      sum += h->GetBinContent(i);
    }
    double per_0pe = sum / total_entries;
    /* find % of 0-PEs */
    TF1 *f0 = new TF1("f0", "[0]*exp(-0.5*((x-[1])/[2])^2)", view_wind[0],
                      view_wind[1]);
    TF1 *f1 = new TF1("f1", "[0]*exp(-0.5*((x-[1])/[2])^2)", view_wind[0],
                      view_wind[1]);
    f0->SetParameter(0, p00);
    f0->SetParameter(1, p01);
    f0->SetParameter(2, p02);
    f1->SetParameter(0, p10);
    f1->SetParameter(1, p11);
    f1->SetParameter(2, p12);
    h->GetXaxis()->SetRangeUser(view_wind[0], view_wind[1]);
    h->Draw();
    h->SetStats(kFALSE);
    f0->Draw("same");
    f1->Draw("same");
    /* add analysis results to plot */
    TPaveText *print_res = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
    ss.str("");
    ss << "# Waveforms " << total_entries;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "Ped. sig.: " << p02;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "P/V: " << pv;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "Res: " << res;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "Adj. SPE peak: " << adj_spe_peak;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "0-PE " << per_0pe << "%";
    print_res->AddText(ss.str().c_str());
    print_res->Draw();
    /* add analysis results to plot */
    c->Update();
    c->Draw();
    c->Print("hist.png");
  }
  else if (fits.size() == 1) {
    TF1 *fit0 = new TF1("fit0", "gaus", fits[0][0], fits[0][1]);
    /* guas function of the form f(x)=p0*exp(-0.5*((x-p1)/p2)^2) */
    h->Fit(fit0, "RQ0");
    double p00 = fit0->GetParameter(0);
    double p01 = fit0->GetParameter(1);
    double p02 = fit0->GetParameter(2);
    double res = p02 / p01;
    TF1 *f0 = new TF1("f0", "[0]*exp(-0.5*((x-[1])/[2])^2)", view_wind[0],
                      view_wind[1]);
    f0->SetParameter(0, p00);
    f0->SetParameter(1, p01);
    f0->SetParameter(2, p02);
    h->GetXaxis()->SetRangeUser(view_wind[0], view_wind[1]);
    h->Draw();
    h->SetStats(kFALSE);
    f0->Draw("same");
    /* add analysis results to plot */
    TPaveText *print_res = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
    ss.str("");
    ss << "# Waveforms " << total_entries;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "Sig.: " << p02;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    ss << "Res: " << res;
    print_res->AddText(ss.str().c_str());
    ss.str("");
    print_res->AddText(ss.str().c_str());
    print_res->Draw();
    /* add analysis results to plot */
    c->Draw();
    c->Update();
    gSystem->ProcessEvents();
  }
}

Results Rooter::fit_spectrum(double wind[2], int view_wind[2],
                          std::vector<std::array<int, 2>> fits, int n_bins) {
  Results results;
  int st = wind[0] / time_res;
  int en = wind[1] / time_res;
  TH1D *h = new TH1D("fit_h", "blank", n_bins, view_wind[0], view_wind[1]);
  int total_entries = tree->GetEntries();
  double integral;
  for (int i = 0; i < total_entries; i++) {
    integral = 0.0;
    tree->GetEntry(i);
    for (int j = st; j < en + 1; j++) {
      integral += (waveform->at(j) - baseline) * polarity;
    }
    h->Fill(integral);
  }
  if (fits.size() == 2) {
    TF1 *fit0 = new TF1("fit0", "gaus", fits[0][0], fits[0][1]);
    TF1 *fit1 = new TF1("fit1", "gaus", fits[1][0], fits[1][1]);
    /* guas function of the form f(x)=p0*exp(-0.5*((x-p1)/p2)^2) */
    h->Fit(fit0, "RQ0");
    h->Fit(fit1, "RQ0+");
    double p00 = fit0->GetParameter(0);
    double p00_err = fit0->GetParError(0);
    double p01 = fit0->GetParameter(1);
    double p01_err = fit0->GetParError(1);
    double p02 = fit0->GetParameter(2);
    double p02_err = fit0->GetParError(2);
    double p10 = fit1->GetParameter(0);
    double p10_err = fit1->GetParError(0);
    double p11 = fit1->GetParameter(1);
    double p11_err = fit1->GetParError(1);
    double p12 = fit1->GetParameter(2);
    double p12_err = fit1->GetParError(2);
    double adj_spe_peak = p11 - p01;
    double adj_spe_peak_err = std::sqrt(p11_err*p11_err + p01_err*p01_err);
    double res = p12 / adj_spe_peak;
    double res_err = std::sqrt((p12_err / p12)*(p12_err / p12) + (adj_spe_peak_err / adj_spe_peak)*(adj_spe_peak_err / adj_spe_peak));
    /* find peak-to-valley */
    h->GetXaxis()->SetRangeUser(p01, p11);
    double min = h->GetMinimum();
    int min_bin = h->GetMinimumBin();
    h->GetXaxis()->SetRange(h->GetMinimumBin(), 5000);
    double max = h->GetMaximum();
    double pv = max / min;
    /* find peak-to-valley */
    /* find % of 0-PEs */
    double sum = 0.0;
    for (int i = 0; i < min_bin + 1; i++) {
      sum += h->GetBinContent(i);
    }
    double per_0pe = sum / total_entries;
    /* find % of 0-PEs */
    std::cout << "# Waveforms " << total_entries << std::endl;
    std::cout << "Ped. sig.: " << p02 << std::endl;
    std::cout << "P/V: " << pv << std::endl;
    std::cout << "Res: " << res << std::endl;
    std::cout << "Adj. SPE peak: " << adj_spe_peak << std::endl;
    std::cout << "0-PE " << per_0pe << "%" << std::endl;
    results.num_waveforms = total_entries;
    results.ped_sig = p02;
    results.ped_sig_err = p02_err;
    results.pv = pv;
    results.res = res;
    results.res_err = res_err;
    results.adj_spe_peak = adj_spe_peak;
    results.adj_spe_peak_err = adj_spe_peak_err;
    results.per_0pe = per_0pe;
  }
  else if (fits.size() == 1) {
    TF1 *fit0 = new TF1("fit0", "gaus", fits[0][0], fits[0][1]);
    /* guas function of the form f(x)=p0*exp(-0.5*((x-p1)/p2)^2) */
    h->Fit(fit0, "RQ0");
    double p00 = fit0->GetParameter(0);
    double p00_err = fit0->GetParError(0);
    double p01 = fit0->GetParameter(1);
    double p01_err = fit0->GetParError(1);
    double p02 = fit0->GetParameter(2);
    double p02_err = fit0->GetParError(2);
    double res = p02 / p01;
    double res_err = std::sqrt((p02_err / p02)*(p02_err / p02) + (p01_err / p01)*(p01_err / p01));
    std::cout << "# Waveforms " << total_entries << std::endl;
    std::cout << "Sig.: " << p02 << std::endl;
    std::cout << "Res: " << res << std::endl;
    results.num_waveforms = total_entries;
    results.ped_sig = p02;
    results.ped_sig_err = p02_err;
    results.res = res;
    results.res_err = res_err;
    results.pv = 0;
    results.adj_spe_peak = 0;
    results.adj_spe_peak_err = 0;
    results.per_0pe = 0;
  }
  return results;
}

void Rooter::view_max_amp(int view_wind[2], int n_bins) {
  TCanvas *c = new TCanvas("max", "Max Amplitudes", 800, 600);
  c->SetLogy();
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Max. Amps.";
  TH1D *h = new TH1D("max_h", ss.str().c_str(), n_bins, view_wind[0], view_wind[1]);
  h->GetYaxis()->SetTitle("Counts");
  h->GetXaxis()->SetTitle("ADC");
  double max_amp = 0;
  double amp = 0;
  int total_entries = tree->GetEntries();
  for (int i = 0; i < total_entries; i++) {
    tree->GetEntry(i);
    max_amp = 0;
    for (int j = 0; j < waveform->size(); j++) {
      amp = (waveform->at(j) - baseline) * polarity;
      if (amp > max_amp)
        max_amp = amp;
    }
    h->Fill(max_amp);
  }
  h->Draw();
  c->Draw();
  c->Update();
  c->Print("max_amp.png");
  gSystem->ProcessEvents();
}

void Rooter::view_waveform(int num) {
  tree->GetEntry(num);
  int waveform_size = waveform->size();
  Double_t x[waveform->size()];
  Double_t y[waveform->size()];
  for (int i = 0; i < waveform_size; i++) {
    x[i] = i * time_res; // convert to ns
    y[i] = (waveform->at(i) - baseline) * polarity;
  }
  TCanvas *c = new TCanvas("full", "Full Waveform", 800, 600);
  TGraph *gr = new TGraph(waveform_size, x, y);
  gr->GetXaxis()->SetTitle("Time [ns]");
  gr->GetXaxis()->SetRangeUser(0, waveform_size * time_res);
  gr->GetYaxis()->SetTitle("ADC");
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Waveform " << num;
  gr->SetTitle(ss.str().c_str());
  gr->Draw();
  c->Update();
  c->Draw();
}

void Rooter::view_waveform_raw(int num) {
  tree->GetEntry(num);
  int waveform_size = waveform->size();
  Double_t x[waveform->size()];
  Double_t y[waveform->size()];
  for (int i = 0; i < waveform_size; i++) {
    x[i] = i * time_res; // convert to ns
    y[i] = waveform->at(i);
  }
  TCanvas *c = new TCanvas("raw", "Raw Waveform", 800, 600);
  TGraph *gr = new TGraph(waveform_size, x, y);
  gr->GetXaxis()->SetTitle("Time [ns]");
  gr->GetXaxis()->SetRangeUser(0, waveform_size * time_res);
  gr->GetYaxis()->SetTitle("ADC");
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Waveform " << num;
  gr->SetTitle(ss.str().c_str());
  gr->Draw();
  c->Update();
  c->Draw();
}

void Rooter::view_waveform_wind(int num, double view_wind[2]) {
  tree->GetEntry(num);
  int waveform_size = waveform->size();
  Double_t x[waveform->size()];
  Double_t y[waveform->size()];
  for (int i = 0; i < waveform_size; i++) {
    x[i] = i * time_res; // convert to ns
    y[i] = waveform->at(i);
  }
  TCanvas *c = new TCanvas("red", "Zoomed Waveform", 800, 600);
  TGraph *gr = new TGraph(waveform_size, x, y);
  gr->GetXaxis()->SetTitle("Time [ns]");
  gr->GetXaxis()->SetRangeUser(view_wind[0], view_wind[1]);
  gr->GetYaxis()->SetTitle("ADC");
  std::stringstream ss;
  ss << "Run " << plot_f_name << " Waveform " << num;
  gr->SetTitle(ss.str().c_str());
  gr->Draw();
  c->Update();
  c->Draw();
}

void Rooter::pre_late_pulsing(std::vector<std::array<int, 2>> winds,
                              int spe_thre, int thre, std::ofstream& ofile,
                              std::string name) {
  int spes = 0;
  int pres = 0;
  int lates = 0;
  int pre_st = winds[0][1] / time_res;
  int pre_en = winds[0][0] / time_res;
  int late_st = winds[1][0] / time_res;
  int late_en = winds[1][1] / time_res;
  int num_entries = tree->GetEntries();
  for (int i = 0; i < num_entries; i++) {
    tree->GetEntry(i);
    for (int j = 0; j < waveform->size(); j++) {
      int adc = (waveform->at(j) - baseline) * polarity;
      if (adc > spe_thre) {
        spes++;
        pres += count_waveform_pulses(thre, j - pre_st, j - pre_en);
        lates += count_waveform_pulses(thre, j + late_st, j + late_en);
        break;
      }
    }
  }
  ofile << name << '\t' << spe_thre << '\t' << static_cast<double>(pres) / spes * 100 << '\t'
        << static_cast<double>(lates) / spes * 100 << '\n';
  // std::cout << "Pre-pulsing: " << static_cast<double>(pres) / spes * 100
  //           << " %\n";
  // std::cout << "Late-pulsing: " << static_cast<double>(lates) / spes * 100
  //           << " %\n";
}

int Rooter::find_average_peak() {
  tree->GetEntry(0);
  int waveform_size = waveform->size();
  std::vector<double> avg_waveform(waveform_size);
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    for (int j = 0; j < waveform_size; j++) {
      avg_waveform[j] += waveform->at(j);
    }
  }
  for (int i = 0; i < waveform_size; i++) {
    avg_waveform[i] /= static_cast<double>(tree->GetEntries());
  }
  double max_amp = 0.0;
  int max_loc;
  for (int j = 0; j < waveform->size(); j++) {
    double amp = (avg_waveform[j] - baseline) * polarity;
    if (amp > max_amp) {
      max_amp = amp;
      max_loc = j;
    }
  }

  return max_loc * time_res;
}
