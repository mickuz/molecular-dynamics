#include "argon.h"

Vector Vector::operator+(const Vector &vector1)
{
    return {this->x + vector1.x, this->y + vector1.y, this->z + vector1.z};
}

Vector Vector::operator-(const Vector &vector1)
{
    return {this->x - vector1.x, this->y - vector1.y, this->z - vector1.z};
}

Vector& Vector::operator+=(const Vector &vector1)
{
    this->x += vector1.x;
    this->y += vector1.y;
    this->z += vector1.z;
    return *this;
}

Vector& Vector::operator-=(const Vector &vector1)
{
    this->x -= vector1.x;
    this->y -= vector1.y;
    this->z -= vector1.z;
    return *this;
}

Vector operator*(const double &n, const Vector &vector1)
{
    return {n*vector1.x, n*vector1.y, n*vector1.z};
}

Vector operator/(const Vector &vector1, const double &n)
{
    return {vector1.x/n, vector1.y/n, vector1.z/n};
}

double Vector::vectorModule()
{
    return sqrt(x * x + y * y + z * z);
}

Parameters::Parameters(char *fileName)
{
    ifstream inputFile(fileName);
	string line;
	getline(inputFile, line);
	n = atoi(line.c_str());
    N = n*n*n;
    getline(inputFile, line);
    m = atoi(line.c_str());
    getline(inputFile, line);
    eps = atoi(line.c_str());
    getline(inputFile, line);
    R = atof(line.c_str());
    getline(inputFile, line);
    f = atof(line.c_str());
    getline(inputFile, line);
    L = atof(line.c_str());
    getline(inputFile, line);
    a = atof(line.c_str());
    getline(inputFile, line);
    T_0 = atoi(line.c_str());
    getline(inputFile, line);
    tau = atof(line.c_str());
    getline(inputFile, line);
    S_o = atoi(line.c_str());
    getline(inputFile, line);
    S_d = atoi(line.c_str());
    getline(inputFile, line);
    S_out = atoi(line.c_str());
    getline(inputFile, line);
    S_xyz = atoi(line.c_str());
    inputFile.close();
}

void State::ComputeInitialPositions(Parameters &par)
{
    Vector b0(par.a, 0, 0);
    Vector b1(par.a/2, par.a*sqrt(3)/2, 0);
    Vector b2(par.a/2, par.a*sqrt(3)/6, par.a*sqrt(2)/sqrt(3));
    Vector* r0 = new Vector[par.N];

    for(int j = 0; j < par.n; j++)
    {
        for (int k = 0; k < par.n; k++)
        {
            for (int l = 0; l < par.n; l++)
            {
                int i = j + k * par.n + l * par.n * par.n;
                r0[i] = (j - (par.n-1) / 2) * b0 + (k - (par.n-1) / 2) * b1 + (l - (par.n-1) / 2) * b2;
            }
        }
    }

    r = r0;
}

void State::ComputeInitialMomenta(Parameters &par)
{
    uniform_real_distribution<double> numberDistribution(0,1);
    uniform_real_distribution<double> signDistribution(-1,1);
    Vector* E_kin = new Vector[par.N];
    Vector* p0 = new Vector[par.N];
    Vector P(0, 0, 0);

    for(int i = 0; i < par.N; i++)
    {
        double xLambda = numberDistribution(generator);
        double yLambda = numberDistribution(generator);
        double zLambda = numberDistribution(generator);
        E_kin[i].setX(-0.5 * k * par.T_0 * log(xLambda));
        E_kin[i].setY(-0.5 * k * par.T_0 * log(yLambda));
        E_kin[i].setZ(-0.5 * k * par.T_0 * log(zLambda));
        int xSign = (signDistribution(generator) < 0) ? -1 : 1;
        int ySign = (signDistribution(generator) < 0) ? -1 : 1;
        int zSign = (signDistribution(generator) < 0) ? -1 : 1;
        p0[i].setX(xSign * sqrt(2 * par.m * E_kin[i].getX()));
        p0[i].setY(ySign * sqrt(2 * par.m * E_kin[i].getY()));
        p0[i].setZ(zSign * sqrt(2 * par.m * E_kin[i].getZ()));
        P += p0[i];
    }

    for(int i = 0; i < par.N; i++)
    {
        p0[i] -= P/par.N;
    }

    p = p0;
}

void State::ComputeForces(Parameters &par)
{
    double distance;
    double lennardjonesPotential = 0;
    double spherePotential = 0;
    double pressure = 0;
    Vector* vanDerWaalsForce = new Vector[par.N];
    Vector* sphereForce = new Vector[par.N];
    Vector* totalForce = new Vector[par.N];

    for(int i = 0; i < par.N; i++)
    {
        distance = r[i].vectorModule();
        if(distance >= par.L)
        {
            spherePotential += 0.5 * par.f * pow(distance - par.L, 2);
            sphereForce[i] = par.f * (par.L - distance) * r[i] / distance;
        }
        pressure += 1/(4*M_PI*par.L*par.L) * sphereForce[i].vectorModule();
        for(int j = i+1; j < par.N; j++)
        {
            distance = (r[i] - r[j]).vectorModule();
            lennardjonesPotential += par.eps * (pow(par.R/distance, 12) - 2 * pow(par.R/distance, 6));
            vanDerWaalsForce[i] += 12 * par.eps * (pow(par.R/distance, 12) - pow(par.R/distance, 6)) * (r[i] - r[j]) / pow(distance, 2);
            vanDerWaalsForce[j] -= 12 * par.eps * (pow(par.R/distance, 12) - pow(par.R/distance, 6)) * (r[i] - r[j]) / pow(distance, 2);
        }
        totalForce[i] = sphereForce[i] + vanDerWaalsForce[i];
    }

    V = lennardjonesPotential + spherePotential;
    F = totalForce;
    P = pressure;
}

void State::UpdateState(Parameters &par)
{
    double temperature = 0;
    double kineticEnergy = 0;

    for(int i = 0; i < par.N; i++)
    {
        p[i] = p[i] + 0.5 * par.tau * F[i];
        r[i] = r[i] + (1.0/par.m) * par.tau * p[i];
    }

    ComputeForces(par);

    for(int i = 0; i < par.N; i++)
    {
        p[i] = p[i] + 0.5 * par.tau * F[i];
        temperature += 1.0 / (3 * par.N * k * par.m) * pow(p[i].vectorModule(), 2);
        kineticEnergy += (0.5/par.m) * pow(p[i].vectorModule(), 2);
    }

    T = temperature;
    H = V + kineticEnergy;
}

void State::SavePositions(Parameters &par, ofstream &outputFile)
{
    for(int i = 0; i < par.N; i++)
    {
        outputFile << "Ar" << i+1 << "	" << r[i].getX() << "	" << r[i].getY() << "	" << r[i].getZ() << endl;
    }
}

void State::SaveData(Parameters &par, ofstream &outputFile, int i)
{
    outputFile << i * par.tau << "  " << H << " " << V << " " << "  " << H-V << "   " << T << " " << P << endl;
}

int main(int argc, char *argv[])
{
    ofstream outputFile, xyzFile;
    Parameters par(argv[1]);
    State argonState;
    argonState.ComputeInitialPositions(par);
    argonState.ComputeInitialMomenta(par);
    argonState.ComputeForces(par);
    outputFile.open(argv[2]);
    outputFile << "t    H   V   K   T   P" << endl;
    xyzFile.open(argv[3]);
    double avgTemperature = 0;
    double avgPressure = 0;
    double avgEnergy = 0;
    for(int s = 0; s < par.getS_o() + par.getS_d(); s++)
    {
        argonState.UpdateState(par);
        if(s % par.getS_out() == 0)
        {
            argonState.SaveData(par, outputFile, s);
        }
        if(s % par.getS_xyz() == 0)
        {
            xyzFile << par.getN() << endl << endl;
            argonState.SavePositions(par, xyzFile);
        }
        if(s >= par.getS_o())
        {
            avgTemperature += argonState.getT();
            avgPressure += argonState.getP();
            avgEnergy += argonState.getH();
        }
    }
    avgTemperature /= par.getS_d();
    avgPressure /= par.getS_d();
    avgEnergy /= par.getS_d();
    cout << endl << "Average temperature: " << avgTemperature << endl << "Average pressure: " << avgPressure << endl << "Average energy: " << avgEnergy;
    outputFile.close();
    xyzFile.close();
    return 0;
}
