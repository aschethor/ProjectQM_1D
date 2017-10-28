
#include "utility.cpp"

using namespace std;

double w = 1;
double m = 1;
double h_bar = 1;
double beta = 10;

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
    vector<vector<double>> MChain;
    vector<double> last = init;
    double p_last = p(last);
    int nBeads = init.size();
    for(int i=0;i<nSteps;i++){
        for(int j=0;j<nBeads;j++){
            vector<double> next = last;
            next[j] += alpha[j]*(-1+2*drand());
            double p_next = p(next);
            if(p_next>p_last||drand()<p_next/p_last){
                p_last = p_next;
                last = next;
            }
            MChain.push_back(last);
        }
    }
    return MChain;
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

/**
 * potential as a function of position x (will be used in probability function p)
 * @param x 
 * @return 
 */
double V(double x){
    //return w*1*(x-1)*(x-1)*(x+1)*(x+1);
    return m*w*w/2*x*x;
}

/**
 * probability to observe system in state with beads at given positions
 * @param positions vector of bead positions
 * @return
 */
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
    
    //set up number of steps, number of beads, initial values, delta
    int nSteps = 10000;
    int N = 10;//number of "beads"
    double delta = 1;
    vector<double> init(N);
    vector<double> v_delta(N);
    for(int i=0;i<N;i++){
        init[i]=0;
        v_delta[i]=delta;
    }
    
    //calculate markov chain
    vector<vector<double>> MCserie = MCMC(&p,init,v_delta,nSteps);
    
    //save markov chain into file
    ofstream myfile;
    myfile.open ("/home/aschethor/Desktop/Uni/11. Semester/Project CQM/doubleWelloutput_qm.csv");
    for(int i=0;i<MCserie.size();){
        double avrg = 0;
        for(int k=0;k<N;k++) {
            for (int j = 0; j < N; j++) {
                avrg += MCserie[i][j];
                myfile<<""<<MCserie[i][j]<<endl;
            }
            i++;
        }
        avrg=avrg/N/N;
        //myfile<<""<<avrg<<endl; //averages classical?!
    }
    myfile.close();
    
    //print out averages
    printf("  N   = %d\n",N);
    printf(" <X>  = %f\n",average(MCserie));
    printf("<X*X> = %f\n",averageOfSquares(MCserie));
    return 0;
}