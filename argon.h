// Created by Michał on 06.11.2018.

#ifndef _argon_h
#define _argon_h
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <ctime>
#include <random>

using namespace std;

static default_random_engine generator{static_cast<unsigned int>(time(nullptr))};
const double k = 8.31E-3;
class State;

class Vector
{
private:
    double x;
    double y;
    double z;

public:
    Vector() : x(0), y(0), z(0) {};
    Vector(double X, double Y, double Z) : x(X), y(Y), z(Z) {};
    Vector operator+(const Vector &vector1);
    Vector operator-(const Vector &vector1);
    Vector& operator+=(const Vector &vector1);
    Vector& operator-=(const Vector &vector1);
    friend Vector operator*(const double &n, const Vector &vector1);
    friend Vector operator/(const Vector &vector1, const double &n);
    void setX(double X) {x = X;}
    void setY(double Y) {y = Y;}
    void setZ(double Z) {z = Z;}
    double getX() {return x;}
    double getY() {return y;}
    double getZ() {return z;}
    double vectorModule();
};

class Parameters
{
friend class State;
private:
    int n; //liczba atomow wzdluz kazdej krawedzi krysztalu
    int N; //calkowita liczba atomow
    int m; //masa atomu
    int eps; //minimum potencjalu
    double R; //odleglosc miedzyatomowa dla ktorej potencjal ma minimum
    double f; //wspolczynnik elastycznosci
    double L; //promien sfery
    double a; //odleglosc miedzy atomami
    int T_0; //zadana temperatura
    double tau; //krok czasowy
    int S_o; //liczba krokow termalizacji
    int S_d; //liczba krokow dynamiki
    int S_out; //czestotliwosc zapisu charakterystyk ukladu do glownego pliku wyjsciowego
    int S_xyz; //czestotliwosc zapisu polozen atomow do pliku XYZ

public:
    Parameters(char *fileName);
    int getN() {return N;}
    int getS_o() {return S_o;}
    int getS_d() {return S_d;}
    int getS_out() {return S_out;}
    int getS_xyz() {return S_xyz;}
};

class State
{
private:
    double V; //energia potencjalna ukladu
    double P; //cisnienie ukladu
    double H; //hamiltonian ukladu
    double T; //temperatura ukladu
    Vector* r; //dynamiczna tablica polozen
    Vector* p; //dynamiczna tablica pedow
    Vector* F; //dynamiczna tablica sil

public:
    void ComputeInitialPositions(Parameters &par);
    void ComputeInitialMomenta(Parameters &par);
    void ComputeForces(Parameters &par);
    void UpdateState(Parameters &par);
    void SavePositions(Parameters &par, ofstream &outputFile);
    void SaveData(Parameters &par, ofstream &outputFile, int i);
    double getT() {return T;}
    double getP() {return P;}
    double getH() {return H;}
};

#endif

