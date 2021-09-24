#include "rooter.hpp"
#include <TApplication.h>
#include <fstream>
#include <iostream>
#include <sstream>

struct Winds {
  int int_wind;
  std::vector<std::array<int, 2>> fit_wind;
};

void write_rate_results();
void write_fit_results(std::ofstream &ofile, std::vector<Rooter *> &runs);
int select_spe_amp(Rooter *r);
Winds select_int_wind_fit_wind(Rooter *r);

int main(int argc, char *argv[]) {
  TApplication app("app", &argc, argv);

  std::vector<Rooter *> runs;
  runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/01.root", 5, 500, 1.44, false, 1400));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/02.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/03.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/04.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/13.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/14.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/15.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/16.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/17.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/18.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/19.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/20.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/21.root", 30, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/24.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/25.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/26.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/30.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/31.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/32.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/36.root", 5, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/38.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/39.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/43.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/45.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/46.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/47.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/51.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/53.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/54.root", 90, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/60.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/62.root", 10, 500));
  runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/65.root", 10, 500, 1.38, true, 1450));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/67.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/69.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/70.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/71.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/74.root", 10, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/81.root", 30, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/82.root", 30, 500));
  // runs.push_back(new Rooter("../../data/summer_2021/ROOT Files/85.root", 10, 500));
  // write_fit_results(file, runs);

  for (int i = 0; i < runs.size(); i++) {
    int view_wind[2] = { -20, 120 };
    runs[i]->view_spectrum(view_wind, 35);
  }

  app.Run();
  return 0;
}

Results analysis(std::vector<Rooter *> &runs) {
  
}

void write_fit_results(std::ofstream &ofile, std::vector<Rooter *> &runs) {
  for (int i = 0; i < runs.size(); i++) {
    int peak = select_spe_amp(runs[i]);
    Winds int_wind = select_int_wind_fit_wind(runs[i]);
    double wind[2];
    int view_wind[] = {-20, 120};
    switch (int_wind.int_wind) {
    case 1:
      wind[0] = 650;
      wind[1] = 680;
      break;
    case 2:
      wind[0] = 640;
      wind[1] = 690;
      break;
    case 3:
      wind[0] = 645;
      wind[1] = 685;
      break;
    }
    Results res = runs[i]->fit_spectrum(wind, view_wind, int_wind.fit_wind, 30);
    std::stringstream ss;
    if (int_wind.fit_wind.size() == 1) {
      ss << int_wind.fit_wind[0][0] << "-" << int_wind.fit_wind[0][1] << ",";
    }
    else {
      ss << int_wind.fit_wind[0][0] << "-" << int_wind.fit_wind[0][1] << ","
         << int_wind.fit_wind[1][0] << "-" << int_wind.fit_wind[1][1];
    }
    ofile << runs[i]->get_red_name() << "," << wind[0] << "-" << wind[1] << ","
          << ss.str() << "," << res.num_waveforms << ","
          << res.ped_sig << "," << res.adj_spe_peak << "," << res.pv << ","
          << res.res << "," << res.per_0pe << "\n";
  }
}

int select_spe_amp(Rooter *r) {
  int view_amp_wind[2];
  if (r->get_amp()) {
    view_amp_wind[0] = 0;
    view_amp_wind[1] = 120;
  }
  else {
    view_amp_wind[0] = 0;
    view_amp_wind[1] = 40;
  }
  int peak = 0;
  r->view_max_amp(view_amp_wind, 40);
  std::cout << "Enter SPE amp. (1/2 of peak): ";
  // std::cout << "Enter SPE amp.: ";
  std::cin >> peak;
  return peak;
}

Winds select_int_wind_fit_wind(Rooter *r) {
  Winds to_return;
  // int int_wind = 0;
  int int_wind = 1;
  int fits;
  int low, high;
  std::array<int, 2> one = {650, 680};
  std::array<int, 2> two = {640, 690};
  std::array<int, 2> three = {645, 685};
  std::vector<std::array<int, 2>> winds;
  winds.push_back(one);
  winds.push_back(two);
  winds.push_back(three);
  int view_spec_wind[] = {-20, 120};
  double spec_wind[] = {650., 680.};
  // r->view_multi_spec(winds, view_spec_wind, 30);
  // r->view_spectrum(spec_wind, view_spec_wind, 30);
  // std::cout << "(1) 650-680\n(2) 640-690\n(3) 645-685\n";
  // std::cout << "Select integration window: ";
  // std::cin >> int_wind;
  to_return.int_wind = int_wind;
  // std::cout << "How many fits to do: ";
  // std::cin >> fits;
  fits = 2;
  for (int i = 0; i < fits; i++) {
    std::cout << "Enter low high: ";
    std::cin >> low >> high;
    std::array<int, 2> fit;
    fit[0] = low;
    fit[1] = high;
    to_return.fit_wind.push_back(fit);
  }
  return to_return;
}

void write_rate_results() {}
