#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>

using namespace std;

bool leerArgumentos(int argc, char *argv[], string& fichero){

	bool f=false;
	string s="";
	for(int i=1; i<argc;i++){
		s = argv[i];
		if(s=="-f"){
			if(i<argc){
				f=true;
				i++;
				fichero=argv[i];
			}
			else{
			cout<<"Debe introducir nombre de fichero despues de -f"<<endl;
			return false;
			}	
		}
		else{
			cout<<"Argumento desconocido: "<<argv[i]<<endl;
			cout<<"Los argumentos que se aceptan son: potter -f ficheroDeEntrada"<<endl;
			return false;
		}
	}
	if(!f){
		cout<<"-f es argumento obligatorio"<<endl;
		return false;
	}
	return true;
}

bool leerFichero(string nombre, int& n, double& T, vector<double>& t, vector<double>& v, vector<unsigned>& m){

	bool abierto=false;
	double aux=0;
	ifstream f(nombre);
	if(f.is_open()){
		abierto=true;
		f>>n;
		f>>T;
		//cout << n << " " << T << endl;
		for(int i=0;i<n;i++){
			f>>aux;			
			t.push_back(aux);
			//cout<<t[i]<<endl;
		}
		for(int i=0;i<n;i++){
			f>>aux;
			v.push_back(aux);
			//cout<<v[i]<<endl;
		}
		for(int i=0;i<n;i++){
			f>>aux;	
			m.push_back(aux);
			//cout<<m[i]<<endl;
		}
	}
	else{
		cout<<"El fichero no existe"<<endl;
		abierto=false;
		f.close();	
	}
	return abierto;
}

double weight(const vector<double>& w, const vector<unsigned>& m, const vector<unsigned>& x) {
	double acc_w = 0.0;
	for (size_t i = 0; i < w.size(); i++)
		acc_w += x[i] * w[i];
	return acc_w;
}

double value(const vector<double>& v, const vector<unsigned>& m, const vector<unsigned>& x) {
	double r = 0.0;
	for (size_t i = 0; i < v.size(); i++) {
		r += v[i] * x[i];//*m[i];
	}
	return r;
}

double tiempo(const vector<double> &t, size_t k, const vector<unsigned> &x) {
	double r = 0.0;
	for (size_t i = 0; i < k; i++) r = t[i] * x[i];
	return r;
}

void knapsack(const vector<double>& v, const vector<unsigned>m, const vector<double>&w, double W,
	size_t k,vector<unsigned>&x, double& best_v, vector<unsigned>& sol) {
	
	if (k == x.size()) {
		if (weight(w, m, x) <= W) {
			double actual_v = value(v, m, x);
			if (actual_v > best_v) {
				best_v = actual_v;
				sol = x;
			}
		}
		return;
	}
	for (unsigned j = 0; j < m[k]; j++) {
		x[k] = j;
		knapsack(v,m, w, W, k + 1, x, best_v,sol);
	}
}

void knapsack(const vector<double>& v,const vector<unsigned>m, const vector<double>& w, double W){
	vector<unsigned>sol(v.size());
	vector<unsigned> x(w.size());
	double best_v = numeric_limits<double>::lowest();
	knapsack(v,m, w, W, 0, x, best_v ,sol);
	cout << "knapsack result: "<<best_v << endl;
	
	cout << "copias result: ";
	for (int i = 0; i < m.size(); i++)
		cout << sol[i]<<" ";
	cout << endl;
}

void imprimir(double resul,double tiempo,const vector<int> copias) {
	cout << resul << endl;
	for (int i = 0; i < copias.size(); i++)
		cout << copias[i]<<" ";
	cout << endl;
	cout << tiempo;
}



void feasiblerec(size_t k,
	const vector<unsigned> m,
	const vector<double>& w, double W,
	vector<unsigned>& x) {
	if (k == x.size()) {
		if (weight(w,m, x) <= W) {
			for (unsigned i = 0; i < x.size(); i++) {
				cout << x[i] << " ";
			}
			cout << endl;
		}
		return;
	}
	for (unsigned j = 0; j < m[k]; j++) {
		x[k] = j;
		feasiblerec(k + 1, m, w, W, x);
	}
}


void feasible(size_t n, 
	const vector<unsigned>& m, 
	const vector<double>& w, double W) {
	vector<unsigned> x(n);
	feasiblerec(0,m,w,W, x);
}




	


//potter -f fichero entrada
int main(int argc, char *argv[]){
	string nombreFichero="";

	//leemos argumentos de entrada
	if(!leerArgumentos(argc,argv, nombreFichero)){
		return 0;
	}
	else{
		int n = -1;
		double T = -1; //num de objetos, maximo de tiempo
		vector<double> v,t;//tiempos
		vector<unsigned> m; // maximo, valor
		
		nombreFichero = "potter_n04r.def";

		if(!leerFichero(nombreFichero,n,T,t,v,m))
			return 0;
		vector<int> copias(n, 0);//contador de copias de cada objeto
		double tiempoTotal = 0;
		double resul = -1;
		//feasible(n,m,t,T);
		knapsack(v,m, t, T);
		//resul = optima(n, T, t, v, m,copias,tiempoTotal);
		imprimir(resul,tiempoTotal,copias);
	}
	return 0;
}
