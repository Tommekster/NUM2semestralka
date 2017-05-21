// # include <iostream>
#include <cmath>
#include <cstdio>
const long double PI = 3.141592653589793238L;

#define OUTPUT_FILE "graph.txt"
#define EXACT_FILE "graph_ex.txt"
//#define EULER

//using namespace std;

class Solution{
  const double a;
  const double h;
  const int num;
  const int dof;
  double *y;
  public:
  Solution(double a, double h, int m,int dof):dof(dof),a(a),h(h),num(m),y(new double [(num+1)*dof]){}
  ~Solution(){delete [] y;}
  double get(int position, int element) {if(position<0) position+=(num+1); return y[dof*position + element];}
  double *get(int position) {if(position<0) position+=(num+1); return &y[dof*position];}
  void set(int position, int element, double value){if(position<0) position+=(num+1); y[dof*position + element] = value;}
  void set(int position, double *value){
    if(position<0) position+=(num+1);
    for(int i =0; i < dof; i++) set(position, i, value[i]);
  }
  void dump(FILE *fout = stdout){
    for(int i = 0; i <= num; i++){
      fprintf(fout, "%10.8f ", i*h+a);
      for(int j = 0; j < dof; j++) fprintf(fout, "%10.8f ", get(i,j));
      fprintf(fout, "\n");
    }
  }
};


class Problem {
  public:
  double gamma1;
  double gamma2;
  double alpha;
  // interval (ia, ib)
  double ia;
  double ib;

  public:
  Problem(double g1, double g2, double a = 0, double b = PI/2)
  :gamma1(g1),gamma2(g2),ia(a),ib(b),alpha((g1+g2)/(b-a)) {}

  int getDegreeOfFreedom(){return 2;}
  void getInitialCondition(double *y0) {
    y0[0] = gamma1;
    y0[1] = alpha;
  }
  void getRightHandSide(double x, double *y, double *y_ /* output */) {
    // -y'' + 16y = 8
    //  y'' = 16*y - 8;
    y_[0] = y[1];         // y' = z;
    y_[1] = 16*y[0] - 8;  // z' = 16*y - 8;
  }
  double aimAt(Solution *solution){
    double F;
    double F_ = 267.7467614837482; // F'(alpha) = y'_alpha (b;alpha) = (exp(2pi)+exp(-2pi))/2

    F = solution->get(-1,1) - gamma2; // F(alpha) = y'(b;alpha) - gamma2
    alpha = alpha - F/F_;  // alpha{k+1} = alpha{k} - F/F'

    return F;
  }

  // exact solution
  double exactY(double x, double c1, double c2){
    return c1*exp(4*x)+c2*exp(-4*x)+0.5;
  }
  double exactY_(double x, double c1, double c2){
    return 4*c1*exp(4*x)-4*c2*exp(-4*x);
  }
  double exactC1(){
    //\frac{\frac{\gamma_2}{4} + \left(\gamma_1 - \frac{1}{2}\right) e^{-2\pi}}{e^{2\pi} + e^{-2\pi}}
    return (gamma2/4.0 - (gamma1-0.5)*exp(-2*PI))/(exp(2*PI)+exp(-2*PI));
    //return exp(-2*PI) * (gamma2/4 - (gamma1-1/2)*exp(-2*PI))/(1+exp(-2*PI));
  }
  double exactC2(double c1){
    return (gamma1 - c1 - 0.5);
  }
  void dumpExactSolution(int num = 100, FILE *fout = stdout){
    double h = (ib-ia)/num;
    double c1 = exactC1();
    double c2 = exactC2(c1);

    for(int i = 0; i <= num; i++){
      double x = i*h+ia;
      fprintf(fout, "%10.8f %10.8f %10.8f ", x, exactY(x,c1,c2), exactY_(x,c1,c2));
      fprintf(fout, "\n");
    }
  }
};

class Solver{
  Problem &problem;
  Solution solution;
  const int num;
  const int dof;
  const double h;
  const double a;
  double * K1;
  double * K2;
  double * K3;
  double * K4;
  double * u;
  double * u_;
  public:
  ~Solver(){
    delete [] K1;
    delete [] K2;
    delete [] K3;
    delete [] K4;
    delete [] u;
    delete [] u_;
  }
  Solver(Problem &p, int m)
    :problem(p)
    ,dof(p.getDegreeOfFreedom())
    ,num(m)
    ,a(p.ia)
    ,h((p.ib-p.ia)/m)
    ,solution(p.ia,(p.ib-p.ia)/m,m,p.getDegreeOfFreedom())
    ,K1(new double [p.getDegreeOfFreedom()])
    ,K2(new double [p.getDegreeOfFreedom()])
    ,K3(new double [p.getDegreeOfFreedom()])
    ,K4(new double [p.getDegreeOfFreedom()])
    ,u(new double [p.getDegreeOfFreedom()])
    ,u_(new double [p.getDegreeOfFreedom()])
    {}
  Solution *getSolution(){return &solution;}
  void solve(){
    problem.getInitialCondition(u);
    solution.set(0,u);
    for(int i = 1; i <= num; i++){
      // K1 = f(x_k,y_k) = f(a+ih, y[])
      problem.getRightHandSide(a + i*h, u,/*output:*/K1);
#ifdef EULER
      for(int j = 0; j < dof; j++) u[j] += h*K1[j];
#else
      // K2 = f(x_k + h/2, y_k + h/2*K1)
      double xk_ = a+i*h + h/2;
      for(int j = 0; j < dof; j++) u_[j] = u[j] + h/2*K1[j];
      problem.getRightHandSide(xk_,u_,/*output:*/K2);
      // K3 = f(x_k + h/2, y_k + h/2*K2)
      for(int j = 0; j < dof; j++) u_[j] = u[j] + h/2*K2[j];
      problem.getRightHandSide(xk_,u_,/*output:*/K3);
      // K4 = f(x_k + h/2, y_k + h/2*K3)
      for(int j = 0; j < dof; j++) u_[j] = u[j] + h/2*K3[j];
      problem.getRightHandSide(xk_,u_,/*output:*/K4);

      for(int j = 0; j < dof; j++) u[j] += h/6*(K1[j] + 2*K2[j] + 2*K3[j] + K4[j]);
#endif
      solution.set(i,u);
    }
  }
};

/*
\[ -y'' + 16 y = 8, \mbox{kde~} x \in (0,\frac{\pi}{2}) \]
\[ y(0) = \alpha\]
\[ y(\frac{\pi}{2}) = \beta \]
*/

int main(int argc, char **argv){
  Problem problem(0.1,0.1);
  Solver solver(problem,1000);

  problem.alpha = -0.4;
  for(int k = 0; k<100; k++){
    printf("alpha(%d) = %.11f\n", k, problem.alpha);
    solver.solve();
    double error = fabs(problem.aimAt(solver.getSolution()));
    if(error < 1e-10) {
      printf("trials: %d\nerror: %.11f\n", k, error);
      break;
    }
  }

  #ifdef OUTPUT_FILE
    FILE *f = fopen(OUTPUT_FILE,"w");
    solver.getSolution()->dump(f);
    fclose(f);
  #else
    printf("Solution:\n");
    solver.getSolution()->dump();
  #endif
  #ifdef EXACT_FILE
    FILE *fe = fopen(EXACT_FILE,"w");
    problem.dumpExactSolution(100,fe);
    fclose(fe);
  #else
    printf("Exact solution:\n");
    problem.dumpExactSolution();
  #endif
}
