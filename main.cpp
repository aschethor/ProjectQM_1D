#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double w = 1;
double m = 10;
double h_bar = 1;
double beta = 10;

double drand(){
    return ((double) rand() / (RAND_MAX));
}

vector<double> add(vector<double>& a,vector<double>& b){
    vector<double> ret(a.size());
    for(int i=0;i<a.size();i++){
        ret[i]=a[i]+b[i];
    }
    return ret;
}

vector<double> multiply(vector<double> a,double b){
    vector<double> ret(a.size());
    for(int i=0;i<a.size();i++){
        ret[i]=a[i]*b;
    }
    return ret;
}

vector<double> multiply(vector<double> a,vector<double> b){
    vector<double> ret(a.size());
    for(int i=0;i<a.size();i++){
        ret[i]=a[i]*b[i];
    }
    return ret;
}

double sum(vector<double> a){
    double ret = 0;
    for(int i=0;i<a.size();i++){
        ret+=a[i];
    }
    return ret;
}

double avrg(vector<double> a){
    return sum(a)/a.size();
}

void print(vector<double> v){
    cout<<"[ ";
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<" ";
    }
    cout<<"]";
}

/**
 *
 * @param p : probability - function
 * @param init : initial values of input parameters
 * @param nSteps : number of marcov chain monte carlo steps
 * @return : vector of size nSteps with values
 */
vector<vector<double >> MCMC(
        double (*p)(vector<double>),
        vector<double> init,
        vector<double> alpha,
        int nSteps){
    vector<vector<double>> ret;
    vector<double> last = init;
    double p_last = p(last);
    int N = init.size();
    for(int i=0;i<nSteps;i++){
        for(int j=0;j<N;j++){
            vector<double> next = last;
            next[j] += alpha[j]*(-1+2*drand());
            double p_next = p(next);
            if(p_next>p_last||drand()<p_next/p_last){
                p_last = p_next;
                last = next;
            }
            ret.push_back(last);
        }
    }
    return ret;
}

double average(vector<vector<double>> values){
    vector<double> sum(values[0].size());
    for(int i=0;i<values.size();i++){
        sum = add(sum,values[i]);
    }
    return avrg(sum)/values.size();
}


double averageOfSquares(vector<vector<double>> values){
    vector<double> sum(values[0].size());
    for(int i=0;i<values.size();i++){
        vector<double> multiple = multiply(values[i],values[i]);
        sum = add(sum,multiple);
    }
    return avrg(sum)/values.size();
}

double V(double x){
    //return w*0.001*(x-1)*(x-1)*(x+1)*(x+1);
    return w*x*x;
}

double p(vector<double> positions){
    double epsilon = beta/positions.size();
    int N = positions.size();
    double rsquares = (positions[N-1]-positions[0])*(positions[N-1]-positions[0]);
    for(int i=1;i<N;i++){
        rsquares+=(positions[i-1]-positions[i])*(positions[i-1]-positions[i]);
    }
    rsquares = m/(2*h_bar*epsilon)*rsquares;
    double potential = 0;
    for(int i=0;i<N;i++){
        potential+=V(positions[i]);
    }
    potential = epsilon*potential;
    double result = exp(-rsquares-potential);
    //printf("P: rsquares: %f potential: %f result: %f avrg-X: %f\n",rsquares,potential,result,avrg(positions));
    return result;
}

int main() {
    srand( (unsigned int)( time(NULL) ) );
    int nSteps = 10000;
    int N = 1;//number of "beeds"
    double alpha = 1;
    vector<double> init(N);
    vector<double> v_alpha(N);
    for(int i=0;i<N;i++){
        init[i]=0;
        v_alpha[i]=alpha;
    }
    vector<vector<double>> MCserie = MCMC(&p,init,v_alpha,nSteps);
    /*for(int i=0;i<MCserie.size();i++){
        print(MCserie[i]);
        printf("\n");
    }*/
    printf("  N   = %d\n",N);
    printf(" <X>  = %f\n",average(MCserie));
    printf("<X*X> = %f\n",averageOfSquares(MCserie));
    return 0;
}