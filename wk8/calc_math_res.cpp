#include <iostream>
#include <string>
#include <fstream>
#include "json.hpp"
#include <cmath>

using namespace std;

void proizv(double x, double v, double w, double g, double& dx_dt, double& dv_dt){
        dx_dt = v;
        dv_dt = -w*w*x - 2*g*v;
}

void stepEiler(double x, double v, double w, double g, double dt, double& x_new, double& v_new){
        double dx_dt;
        double dv_dt;
        proizv(x, v, w, g, dx_dt, dv_dt);
        x_new = x + dt*dx_dt;
        v_new = v + dt*dv_dt;
}


void stepHoin(double x, double v, double w, double g, double dt, double& x_new, double& v_new){
        double dx_dt1;
        double dv_dt1;
        proizv(x, v, w, g, dx_dt1, dv_dt1);
        double x_pred = x + dt*dx_dt1;
        double v_pred = v + dt*dv_dt1;

        double dx_dt2, dv_dt2;
        proizv(x_pred, v_pred, w, g, dx_dt2, dv_dt2);
        x_new = x + dt*(dx_dt1 + dx_dt2)/2.0;
        v_new = v + dt*(dv_dt1 + dv_dt2)/2.0;
}

void steprk(double x, double v, double w, double g, double dt, double& x_new, double& v_new){
        double dx_dt1, dv_dt1;
        proizv(x, v, w, g, dx_dt1, dv_dt1);

        double dx_dt2, dv_dt2;
        proizv(x + 0.5*dx_dt1*dt, v + 0.5*dv_dt1*dt, w, g, dx_dt2, dv_dt2);

        double dx_dt3, dv_dt3;
        proizv(x + 0.5*dx_dt2*dt, v + 0.5*dv_dt2*dt, w, g, dx_dt3, dv_dt3);

        double dx_dt4, dv_dt4;
        proizv(x + dx_dt3*dt, v + dv_dt3*dt, w, g, dx_dt4, dv_dt4);

        x_new = x + dt*(dx_dt1 + 2*dx_dt2 + 2*dx_dt3 + dx_dt4)/6.0;
        v_new = v + dt*(dv_dt1 + 2*dv_dt2 + 2*dv_dt3 + dv_dt4)/6.0;
}

int main(int argc, char** argv){
        ifstream in(argv[1]);
        nlohmann::json j;
        in >> j;

        double x = j.at("x0").get<double>();
        double v = j.at("v0").get<double>();
        double w = j.at("w").get<double>();
        double g = j.at("g").get<double>();
        double dt = j.at("dt").get<double>();
        double t_end = j.at("t_end").get<double>();
        string method = j.at("method").get<string>();
        string name_out = j.at("name_out").get<string>();
        double x_new = 0.0;
        double v_new = 0.0;

        ofstream fout(name_out, ios::binary | ios::trunc);

        for(double t = 0; t < t_end; t += dt){
                fout.write(reinterpret_cast<const char*>(&x), sizeof(x));
                fout.write(reinterpret_cast<const char*>(&v), sizeof(v));
                fout.write(reinterpret_cast<const char*>(&t), sizeof(t));
                if(method == "eiler"){
                        stepEiler(x, v, w, g, dt, x_new, v_new);
                        x = x_new;
                        v = v_new;
                }else if(method == "hoin"){
                        stepHoin(x, v, w, g, dt, x_new, v_new);
                        x = x_new;
                        v = v_new;
                }else if(method == "rk"){
                        steprk(x, v, w, g, dt, x_new, v_new);
                        x = x_new;
                        v = v_new;
                }
        }
        fout.close();
}
