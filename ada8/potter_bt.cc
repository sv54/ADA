#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<numeric>

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

double voraz(int n, double T, const vector<double>& t, const vector<double>& v, const vector<unsigned>& m){//, double& tiempo) {
	vector <double> vt(n);

	for (int i = 0; i < n; i++) {
		vt[i] = i;//
	}

	sort(vt.begin(), vt.end(), [&v, &t](size_t x, size_t y) {
		return v[x] / t[x] > v[y] / t[y]; });

	double total = 0;
	for (auto i : vt) {
		for (int j = m[i]; j > 0; j--) {

			if (t[i] * j <= T) {
				total += j * v[i];
				//copias[i] = copias[i] + j;
				//tiempo += t[i] * j;
				T -= t[i] * j;
				j = j - j;
				i++;
			}
		}
	}

	return total;
}

ostream& operator<<(ostream& os, const vector<unsigned>& vect) {
	for (unsigned i = 0; i < vect.size(); i++) {
		os << vect[i] << " ";
	}
	return os;
}

double weight(const vector<double>& w, size_t k, const vector<unsigned>& x) {
	double acc_w = 0.0;
	for (size_t i = 0; i < k; i++)
		acc_w += x[i] * w[i];
	return acc_w;
}

double value(const vector<double>& v, const vector<unsigned>& x) {
	double r = 0.0;
	for (size_t i = 0; i < v.size(); i++) {
		r += v[i] * x[i];
	}
	return r;
}

vector<double> knapsack_W(
	const vector<double>& v, // values
	const vector<double>& w, // weights
	double W // knapsack weight limit
) {
	vector<unsigned> idx(w.size());
	for (unsigned i = 0; i < idx.size(); i++) idx[i] = i;
	sort(idx.begin(), idx.end(), [v, w](unsigned x, unsigned y) {
		return v[x] / w[x] > v[y] / w[y]; });

	vector<double> x(w.size());
	double acc_v = 0.0;
	for (unsigned i = 0; i < idx.size(); i++) {
		if (w[idx[i]] > W) {
			acc_v += W / w[idx[i]] * v[idx[i]];
			x[idx[i]] = W / w[idx[i]];
			W = 0.0;

		}
		else {
			acc_v += v[idx[i]];
			W -= w[idx[i]];
			x[idx[i]] = 1.0;

		}

	}
	return x;

}


double knapsack_c(const vector<double>& v, // values
				const vector<double>& w, // weights
				double W) {
	vector<unsigned> idx(w.size()); // objects sorted by value density
	for (unsigned i = 0; i < idx.size(); i++) idx[i] = i;

	sort(idx.begin(), idx.end(),
		[v, w](unsigned x, unsigned y) {
			return v[x] / w[x] > v[y] / w[y];
		}
	);
	double acc_v = 0.0;
	for (unsigned i = 0; i < idx.size(); i++) {
		if (w[idx[i]] > W) {
			acc_v += W / w[idx[i]] * v[idx[i]];
			break;

		}
		acc_v += v[idx[i]];
		W -= w[idx[i]];

	}
	return acc_v;
}


double add_rest(const vector<double>& v, const vector<unsigned>& m,size_t k) {
	double res = 0.0;
	for (size_t i = k; i < v.size(); i++)
		res += v[i]*m[i];
	return res;
}

void knapsack(const vector<double>& v, const vector<unsigned>m, const vector<double>& w, double W,
	size_t k, vector<unsigned>& x, double& best_v, double& acc_w, double& acc_v) {


	if (k == x.size()) {
		best_v = max(best_v, acc_v);
		/*
		if (weight(w, k, x) <= W) {
			double actual_v = value(v, x);
			if (actual_v > best_v) {
				best_v = actual_v;
				sol = x;
			}
		}*/ //si queremos saber la solucion
		return;
	}
	
	for (unsigned j = 0; j < m[k]; j++) {
		x[k] = j;
		double present_w = acc_w + x[k] * w[k];
		double present_v = acc_v + x[k] * v[k];
		if (present_w<= W &&
			present_v + knapsack_c(v,w,W-present_w)>best_v)
			knapsack(v,m, w, W, k + 1, x, best_v, present_w, present_v);
	}
	
}

void knapsack(const vector<double>& v,const vector<unsigned>m, const vector<double>& w, double W){
	vector<unsigned>sol(v.size());
	vector<unsigned> x(w.size());

	vector<size_t> idx(v.size()); // index vector
	iota(begin(idx), end(idx), 0);
	sort(begin(idx), end(idx),
		[&v, &w](size_t i, size_t j) {
			return v[i] / w[i] > v[j] / w[j];
		}
	);
	vector<double> s_v(v.size()), s_w(w.size());

	for (size_t i = 0; i < v.size(); i++) {
		s_v[i] = v[idx[i]]; // sorted values
		s_w[i] = w[idx[i]]; // sorted weights

	}


	double auxCeroDouble = 0;
	double best_v = voraz(v.size(),W,s_w,s_v,m);
	cout << "voraz result: "<<best_v << endl;
	knapsack(s_v,m, s_w, W, 0, x, best_v, auxCeroDouble, auxCeroDouble);
	cout << "knapsack result: "<<best_v << endl;
	
	cout << "copias result: " << sol;
	/*for (int i = 0; i < m.size(); i++)
		cout << sol[i]<<" ";*/
	cout << endl;
}

void imprimir(double resul,double tiempo,const vector<int> copias) {
	cout << resul << endl;
	for (int i = 0; i < copias.size(); i++)
		cout << copias[i]<<" ";
	cout << endl;
	cout << tiempo;
}

/*
void feasiblerec(size_t k,
	const vector<unsigned> m,
	const vector<double>& w, double W,
	vector<unsigned>& x) {
	if (k == x.size()) {
		if (weight(w, x) <= W) {
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


*/

	
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
		

		nombreFichero = "potter_n15.def";

		if(!leerFichero(nombreFichero,n,T,t,v,m))
			return 0;
		vector<int> copias(n, 0);//contador de copias de cada objeto
		double tiempoTotal = 0;
		double resul = -1;
		//feasible(n,m,t,T);
		clock_t start = clock();
		knapsack(v,m, t, T);
		clock_t end = clock();
		cout << (double(end - start) / ((clock_t)1000)) << "s" << endl;
		//resul = optima(n, T, t, v, m,copias,tiempoTotal);
		imprimir(resul,tiempoTotal,copias);
	}
	return 0;
}
