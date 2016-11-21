// Copyright (c) 2016 Richard Glennie, University of St Andrews 
// 
// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files, to deal in the
// software without restriction, including without limitation the right to use,
// copy, modify, publish, distribute, sublicense, and/or sell copies of the software,
// and to permit persons to whom the software is furnished to do so, subject to the 
// following conditions: 
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the software. 
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS  OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE THE USE OR OTHER DEALINGS
// IN THE SOFTWARE
//
// ============================================================================
//
// Project: Distance Sampling with movement simulation engine 
// File contains cpp code to be included using Rcpp
//
// ============================================================================
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <fstream>
#include <cmath> 
#include <armadillo> 
#include <errno.h> 
#include <RcppArmadillo.h> 
#include <vector> 
#include <string> 
#include <sstream> 

using arma::mat; 
using arma::vec;
using arma::rowvec;
using arma::field; 
using arma::cx_mat;
using arma::cx_cube; 
using arma::randu;
using arma::randn; 

using std::cout;
using std::string; 

// Class Declarations
// ============================================================================
class SurveyRegion;
class LineTransect;
class Animal;
class Observer; 

// The survey region is the area the population lives in and contains 
// the total population. 
// An region is specified by its width and length (it is assumed to be rectangular). 
class SurveyRegion { 
  public:
    // constructors
    SurveyRegion(); 
    SurveyRegion(double width, double length) : length_(length), width_(width) {} 
    SurveyRegion(vec size) : length_(size(1)), width_(size(0)){}

    // accessors 
    const double width() const {return width_;}
    const double length() const {return length_;} 

  private: 
    double length_; 
    double width_; 
};

// Line transect is tranversed by the observer.
// It is specified by its length and its width (full width, not half-width). 
// For the simulation, it is assumed the transect is placed (w / 2, 0) where 
// w is the width of the survey region. As animals are on average uniformly
// placed with respect to the line, where the line does not matter. 
// Adding a starting coordinate for the line transect is an obvious extension.  
class LineTransect {
  public:  
    // Default constructor is empty (does nothing)
    LineTransect(); 
    LineTransect(double width, double length) : length_(length), width_(width) {} 
    LineTransect(vec size) : length_(size(1)), width_(size(0)) {} 

    // Accessor functions 
    const double length() const {return length_;} 
    const double width() const {return width_;}

  private: 
    double length_; 
    double width_; 
};

// Animals are assumed to move according to behaviour: 
//   0 = no movement
//   1 = random direction, constant speed 
//   2 = OU process with home center
// Animals also have a accumulated hazard and a threshold hazard. An animal
// accumulates hazard by being surveyed by an observer. When the accumulated 
// hazard exceeds the threshold, the animal is seen by an observer. 
// The accumulation is setup for only one observer. 
// Animals can only accumulate hazard when they are available for detection. 
// It is assumed that animals are reflected off the boundary of the survey region. 
class Animal {
  public:
    // constructors
    Animal(){}  
    // place animal at (x,y), not moving
    Animal(double x, double y) : behaviour_(0), available_(true) {
      x_ = x; 
      y_ = y;
      threshold_hazard_ = -log(arma::as_scalar(randu(1))); 
      accumulated_hazard_ = 0;  
    }
    // place animal at (x,y) moving in direction at velocity v
    Animal(double x, double y, double v, double direction) : x_(x), y_(y), v_(v),
    direction_(direction), behaviour_(1), available_(true) {
      threshold_hazard_ = -log(arma::as_scalar(randu(1)));
      accumulated_hazard_ = 0;  
    }
    // place animal at home center (x_home, y_home) and 
    // move as OU process with parameters theta = (tau, c)
    // See Gillespie et al. (1996) for parameterisation.
    Animal(double x_home, double y_home, vec theta) : x_(x_home), y_(y_home), behaviour_(2), available_(true), x_home_(x_home), y_home_(y_home), movement_parameter_(theta) {
      threshold_hazard_ = -log(arma::as_scalar(randu(1))); 
      accumulated_hazard_ = 0; 
    } 
    // accessor functions 
    const double x() const {return x_;}
    const double y() const {return y_;} 
    const double v() const {return v_;} 
    const double direction() const {return direction_;} 
    const int behaviour() const {return behaviour_;} 
    const double x_home() const {return x_home_;} 
    const double y_home() const {return y_home_;} 
    const bool available() const {return available_;}
    const double hazard() const {return accumulated_hazard_;}
    const double threshold_hazard() const {return threshold_hazard_;} 
    const vec movement_parameter() const {return movement_parameter_;}  

    // mutator functions 
    void set_x(double x){x_ = x;} 
    void set_y(double y){y_ = y;} 
    void set_v(double v){v_ = v;} 
    void set_direction(double direction){direction_ = direction;}
    void set_behaviour(int behaviour){behaviour_ = behaviour;} 
    void set_home(double x_home, double y_home){x_home_ = x_home; y_home_ = y_home;}  
    void set_available(bool available){available_ = available;} 
    void set_movement_parameter(vec movement_parameter){movement_parameter_ = movement_parameter;}  
    void reset_threshold() {threshold_hazard_ = -log(arma::as_scalar(randu(1)));} 
    void reset_hazard(){accumulated_hazard_ = 0;}
    void accumulate_hazard(const Observer observer, double dt);  
    void move(const SurveyRegion& region, double dt);  
  private: 
    double x_;
    double y_; 
    double v_;
    double direction_; 
    int behaviour_; 
    double x_home_; 
    double y_home_; 
    bool available_;
    double accumulated_hazard_;
    double threshold_hazard_;  
    vec movement_parameter_;  
}; 

// Observer moves along transect, detecting animals.  
// Assumes observer starts at beginning of transect
// Detection functions is 2D hazard, so detection_parameter = (c, b) 
// See advanced distance sampling book p350 by Buckland et al.   
class Observer {
  public:
    Observer(); 
    Observer(double v, vec detection_parameter) : x_(0), y_(0), v_(v), detection_parameter_(detection_parameter){}   

    // accessor functions 
    const double x() const {return x_;}
    const double y() const {return y_;} 
    const double speed() const {return v_;} 
    const vec detection_parameter() const {return detection_parameter_;} 

    // mutator functions
    void set_x(double x){x_ = x;} 
    void set_y(double y){y_ = y;} 
    void set_speed(double v){v_ = v;}
    void set_detection_parameter(vec theta) {detection_parameter_ = theta;} 

    void move(double dt); 
    bool detect(Animal& animal, double dt); 

  private: 
    double x_;
    double y_; 
    double v_;
    vec detection_parameter_;  
}; 

// SurveyDat is a data structure used to package information
// together for brevity.
struct SurveyDat {
  int population_size; 
  vec region_size;
  vec line_size; 
  int behaviour; 
  vec movement_parameter;   
  double observer_speed;
  vec detection_parameter;  
  int num_transects; 
};



// Class Member function definitions 
// ============================================================================

// Accumulates hazard of detection based on relative distance between 
// observer and animal. 
// Inputs: 
//   observer: an Observer object, the observer who is surveying 
//   dt: time step for simulation 
void Animal::accumulate_hazard(const Observer observer, double dt) {
  double rel_x = x() - observer.x(); 
  double rel_y = y() - observer.y(); 
  // animals are assumed to be detectable iff they are available and 
  // are positioned in FRONT of the observer
  // (this was a somewhat arbitrary choice, no reason that y > 0 is needed) 
  if (available() & rel_y > 0) {
    double r0 = rel_x * rel_x + rel_y * rel_y;
    double y1 = rel_y - observer.speed() * dt;
    // if y1 < 0 but y > 0, animal lies in front of observer now 
    // but observeer will pass the animal abeam in this time step, so 
    // only survey up until y1 = 0  
    if (y1 < 0) y1 = 0; 	
    double r1 = rel_x * rel_x + y1 * y1;
    // animals at same position as observer causes overflow in hazard, 
    // so handle seperately.  
    if (r1 < 1e-10) accumulated_hazard_ = arma::datum::inf; 
    else {  
      double alpha = observer.detection_parameter()(0); 
      double beta = observer.detection_parameter()(1); 
      double abeta = 0.5 * (beta - 1); 
      // hazard is integral over time, if x is small or 
      // abeta is small, then integral differs from usual case
      if (fabs(rel_x) < 1e-10) {
        if (fabs(abeta) < 1e-10) accumulated_hazard_ += alpha * 0.5 * (log(r1) - log(r0)); 
        else accumulated_hazard_ += beta / (alpha - 1) * (1 / pow(r1, abeta) - 1 / pow(r0,abeta)); 
      }
      else {
        double hazard = R::pbeta(rel_x * rel_x / r1, abeta, 0.5, 1, 0) - R::pbeta(rel_x * rel_x / r0, abeta, 0.5, 1, 0); 
        hazard *= R::beta(abeta, 0.5) * alpha / (2 * pow(fabs(rel_x), beta - 1)); 
        accumulated_hazard_ += hazard; 
      }
    }      
  }
}

// Move animal over a time interval of length dt
// Inputs: 
//   region: a SurveyRegion object
//   dt: time step 
void Animal::move(const SurveyRegion& region, double dt) {
  double new_x = x(); 
  double new_y = y();
  double new_direction = direction();  
  switch (behaviour()) { 
    // straight-line, random-direction, constant-speed movement
    case 1: 
      new_x = x() + v() * cos(direction()) * dt; 
      new_y = y() + v() * sin(direction()) * dt; 
      break; 
      // OU home-range process 
    case 2:
      double var = movement_parameter()(1) * movement_parameter()(0) / 2; 
      var *= 1 - exp(-2 * dt / movement_parameter()(0)); 
      var = sqrt(var); 
      new_x = x_home() + (x() - x_home()) * exp(-dt / movement_parameter()(0)); 
      new_x += var * as_scalar(randn(1)); 
      new_y = y_home() + (y() - y_home()) * exp(-dt / movement_parameter()(0)); 
      new_y += var * as_scalar(randn(1)); 
      break; 
  }
  // check if animal has moved outside region, if so, reflect off boundary 
  // and reverse direction 
  if (new_x > region.width()) {
    new_x = 2 * region.width() - new_x;
    if (behaviour() == 1) new_direction = fmod(new_direction + M_PI, 2 * M_PI);  
  }
  if (new_y > region.length()) {
    new_y = 2 * region.length() - new_y; 
    if (behaviour() == 1) new_direction = fmod(new_direction + M_PI, 2 * M_PI);  
  }
  if (new_x < 0) {
    new_x *= -1;
    if (behaviour() == 1) new_direction = fmod(new_direction + M_PI, 2 * M_PI);  
  }
  if (new_y < 0) {
    new_y *= -1; 
    if (behaviour() == 1) new_direction = fmod(new_direction + M_PI, 2 * M_PI);  
  }
  set_x(new_x);
  set_y(new_y);
  if (behaviour() == 1) set_direction(new_direction); 
} 

void Observer::move(double dt) {
  set_y(y() + speed() * dt); 
}

bool Observer::detect(Animal& animal, double dt) {
  // if hazard > threshold already, the animal has already been detected, 
  // you don't detect it again. 
  if (animal.hazard() > animal.threshold_hazard()) return false;
  animal.accumulate_hazard(*this,  dt); 
  // if hazard > threshold now, animal is seen at this time
  if (animal.hazard() > animal.threshold_hazard()) return true; 
  return false; 	
}  

// Standalone Functions 
// ============================================================================

// Respawn animals for a new transect
// Inputs: 
//   population: pointer to array of animals
//   num_animals: number of animals in population
//   region: SurveyRegion object 
void SpawnAnimals(std::vector<Animal>& population, const int& num_animals, const SurveyRegion& region) {
  // animals are assumed to be randomly distributed about the region
  // (and hence line)  
  vec randx = randu(num_animals) * region.width(); 
  vec randy = randu(num_animals) * region.length(); 
  for (int animal = 0; animal < num_animals; ++animal) {
    population[animal].set_x(randx(animal));
    population[animal].set_y(randy(animal)); 
    population[animal].reset_hazard(); 
    population[animal].reset_threshold();
    population[animal].set_available(true);
    switch (population[animal].behaviour()) {
      case 1:
        population[animal].set_direction(as_scalar(randu(1)) * M_PI * 2);
        population[animal].set_v(population[animal].movement_parameter()(0)); 
        break; 
      case 2: 
        // home ranges are assumed to be randomly distributed
        population[animal].set_home(randx(animal), randy(animal));	
        break;
    }    
  }
}

// Respawn observer for new transect. Observer is placed at (width/2, 0) 
//  coordinates.  
// Inputs: 
//   observer: Observer object to be respawned 
//   region: SurveyRegion object 
void SpawnObserver(Observer& observer, const SurveyRegion& region) {
  observer.set_x(region.width() / 2); 
  observer.set_y(0); 
}	

// Simulate distance sampling survey with moving animals
// Inputs: 
//   simdat: SimulationDat object
//   dt: time step 
// Output: 
//   file "data.txt" is written to mcds directory with tabulated fields
//   "transect", "effort", "perpendicular distance" for each recorded detection.
void SimulateSurvey(SurveyDat survey_dat, double dt) {
  // create objects involved
  SurveyRegion region(survey_dat.region_size); 
  LineTransect line(survey_dat.line_size);
  int num_animals = survey_dat.population_size;  
  std::vector<Animal> population(num_animals);
  Observer observer(survey_dat.observer_speed, survey_dat.detection_parameter); 
  for (int animal = 0; animal < num_animals; ++animal) {
    population[animal].set_behaviour(survey_dat.behaviour);
    population[animal].set_movement_parameter(survey_dat.movement_parameter); 
  }

  // open file to write data to
  std::ofstream dat("./mcds/data.txt");

  // perform survey 
  double num_timesteps = floor(line.length() / (dt * observer.speed()));
  vec randx(num_animals); 
  vec randy(num_animals);
  int num_detected;  
  for (int transect = 0; transect < survey_dat.num_transects; ++transect) {
    SpawnAnimals(population, num_animals, region); 
    SpawnObserver(observer, region);  
    num_detected = 0; 
    for (int t = 0; t < num_timesteps; ++t) {
      for (int animal = 0; animal < num_animals; ++animal) { 
        // animal detection is recorded if animal is seen INSIDE transect only
        // (but notice this does not ignore animals seen outside transect, it merely doesn't record them)  
        if (observer.detect(population[animal], dt) && fabs(population[animal].x() - observer.x()) <= line.width()/2.0) {
          dat << transect << "\t"  << line.length() << "\t" << fabs(population[animal].x() - observer.x()) << "\n";
          ++num_detected; 
        } 
        population[animal].move(region, dt); 
      } 
      observer.move(dt); 
    } 
    // if no animals seen on a transect, write a line to indicate this for MCDS 
    if (num_detected == 0) dat << transect << "\t" << line.length() << "\n"; 
  }
}

// Fit distance sampling model
// Assumes there is a file named "data.csv" with fields 
// transect, perpendicular distance.
// This function runs the external MCDS engine; the engine must 
// be in the directory "./mcds", and the command file for the 
// MCDS engine must be correctly filled out. 
// Current runs MCDS using WINE! 
// Outputs:
//   writes a file named "output.txt" with the model results
void FitModel() {
  int r = system("cd mcds; wine MCDS 0, command.txt"); 
}


// Save the data and fitted model of simulated survey to 
// a subdirectory. 
// The function renames "data.txt" to "../results/data_<num>.csv"
// and output from mcds "output.txt" to "../results/output_<num>.txt"
// where "<num>" is the simulation number. 
// Note this function assumes the directiory "results" exists.
// Input: 
//   num: the number to appended to each file name 
//   save_log: if TRUE, the MCDS log files are saved also (to check for convergence problems)
void Clean(int num, bool save_log = false) {
  std::stringstream cmd; 
  cmd << "mv ./mcds/data.txt ../results/data_" << num << ".txt"; 
  int r = system(cmd.str().c_str()); 
  cmd.str(""); 
  cmd << "mv ./mcds/output.txt ../results/output_" << num << ".txt"; 
  r = system(cmd.str().c_str());  
  if (save_log) {
    cmd.str("");
    cmd << "mv ./mcds/log.txt ../results/log_" << num << ".txt";
    r = system(cmd.str().c_str());  
  }
}

// Perform Simulation study
// Inputs:
//   num_simulations: number of simulated surveys to do
//   parameter: vector of detection parameters (c,d) and 
//     if behaviour=1 (speed, angle sd)
//     if behaviour=2 (tau, sigma) of OU process 
//   simulation_dat: vector containing (in this order)  
//     population_size: size of population to simulate
//     region_size: (width, length) of rectangular study region
//     line_size: (width, length) of a single transect 
//     behaviour: 0(no movement), 1(straight-line movement), 2(home-range movement) 
//     observer_speed: speed of the observer 
//     num_transects: number of transects per surveyed simulation 
//   dt: time step
//   seed: random seed (set in Armadillo) 
// [[Rcpp::export]] 
void Simulate(int num_simulations, vec parameter, vec simulation_dat, double dt = 1, int seed = -1) {
  if (seed > 0) arma::arma_rng::set_seed(seed); 
  // unpack simulation data from input and then
  // package into a data structure  
  int population_size = simulation_dat(0); 
  vec region_size = simulation_dat.rows(1,2); 
  vec line_size = simulation_dat.rows(3,4); 
  int behaviour = simulation_dat(5); 
  double observer_speed = simulation_dat(6); 
  int num_transects = simulation_dat(7);
  vec detection_parameter = parameter.rows(0,1); 
  vec movement_parameter = arma::zeros<vec>(2); 
  switch (behaviour) {
    case 1: 
      movement_parameter(0) = parameter(2); 
      break; 
    case 2:
      movement_parameter = parameter.rows(2,3);
      break;  
  }
  SurveyDat survey_dat = {population_size, region_size, line_size, behaviour, movement_parameter, observer_speed, detection_parameter, num_transects}; 	

  for (int sim = 0; sim < num_simulations; ++sim) {
    cout << sim + 1 << " / " << num_simulations << "\n"; 
    SimulateSurvey(survey_dat, dt);
    FitModel(); 
    Clean(sim + 1); 	
  }
}	
