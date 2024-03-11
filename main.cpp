//
// Created by Vitto Resnick on 2/21/24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

// Check if a number is an integer
bool is_integer(double k){
    return floor(k) == k;
}

// Factorial function
int factorial(int n){
    int res = 1,i;
    for (i=1;i<=n;i++){
        res *= i;
    }
    return res;
}

// Double factorial
int dfact(int n){
    int i;double res=1.0;
    for(i=n;i>=1;i-=2){
        res *=i;
    }
    return res;
}

// Binomial Coefficient Example: m choose n
int binomCoeff(int m, int n){
    return factorial(m)/(factorial(n)*factorial(m-n));
}

// Calculation one of three directional components for SAB (overlap integral of prod of 2 gauss)
double SxABComp(double XA, double XB, double alpha, double beta, double lA, double lB){

    double P  = exp(-alpha*beta*pow(XA-XB,2)/(alpha + beta)); // Calculate prefactor
    double XP = (alpha*XA + beta*XB)/(alpha + beta);
    double doubleSum = 0; // Initialize double sum

    // Compute  double sound
    for     (int i = 0; i < lA+1; i++){
        double innerSum = 0;
        for (int j = 0; j < lB+1; j++){
            if ((i+j)% 2 == 0){ // Only do even i+j terms
                double summand = binomCoeff(lA,i)*binomCoeff(lB,j)*dfact(i+j-1)*pow(XP-XA,lA-i)*pow(XP-XB,lB-j)/pow(2*(alpha+beta),(i+j)/2);
                innerSum += summand;
            }
        }
        doubleSum += innerSum;
    }
    return P*sqrt(M_PI/(alpha + beta))*doubleSum;
}

// Find normalization constant for a given k primitive gauss
double NormConst(double X, double Y, double Z, double alpha_k, double l, double m, double n){
    double SAA = 1; // Initialize SAA product
    SAA *= SxABComp(X,X,alpha_k,alpha_k,l,l); // Compute SxAA
    SAA *= SxABComp(Y,Y,alpha_k,alpha_k,m,m); // Compute SyAA
    SAA *= SxABComp(Z,Z,alpha_k,alpha_k,n,n); // Compute SzAA

    // normalization constants are defined such that the overlap of each primitive with itself is equal to one
    double N_k_lmn = 1/sqrt(SAA);
    return N_k_lmn;
}

// Process basis data from basis files
vector<vec> ProcessBasisData(string H_STO3G, string C_STO3G, string N_STO3G, string O_STO3G, string F_STO3G){
    double H_alpha_1s_1,H_alpha_1s_2,H_alpha_1s_3;
    double H_d_1s_1,H_d_1s_2,H_d_1s_3;

    double C_alpha_2s_2p_1,C_alpha_2s_2p_2,C_alpha_2s_2p_3;
    double C_d_2s_1       ,C_d_2s_2       ,C_d_2s_3       ;
    double C_d_2p_1       ,C_d_2p_2       ,C_d_2p_3       ;

    double N_alpha_2s_2p_1,N_alpha_2s_2p_2,N_alpha_2s_2p_3;
    double N_d_2s_1       ,N_d_2s_2       ,N_d_2s_3       ;
    double N_d_2p_1       ,N_d_2p_2       ,N_d_2p_3       ;

    double O_alpha_2s_2p_1,O_alpha_2s_2p_2,O_alpha_2s_2p_3;
    double O_d_2s_1       ,O_d_2s_2       ,O_d_2s_3       ;
    double O_d_2p_1       ,O_d_2p_2       ,O_d_2p_3       ;

    double F_alpha_2s_2p_1,F_alpha_2s_2p_2,F_alpha_2s_2p_3;
    double F_d_2s_1       ,F_d_2s_2       ,F_d_2s_3       ;
    double F_d_2p_1       ,F_d_2p_2       ,F_d_2p_3       ;


    ifstream inputFile0(H_STO3G);
    inputFile0 >> H_alpha_1s_1     >> H_d_1s_1;
    inputFile0 >> H_alpha_1s_2     >> H_d_1s_2;
    inputFile0 >> H_alpha_1s_3     >> H_d_1s_3;
    inputFile0.close();

    ifstream inputFile1(C_STO3G);
    inputFile1 >> C_alpha_2s_2p_1 >> C_d_2s_1 >> C_d_2p_1;
    inputFile1 >> C_alpha_2s_2p_2 >> C_d_2s_2 >> C_d_2p_2;
    inputFile1 >> C_alpha_2s_2p_3 >> C_d_2s_3 >> C_d_2p_3;
    inputFile1.close();

    ifstream inputFile2(N_STO3G);
    inputFile2 >> N_alpha_2s_2p_1 >> N_d_2s_1 >> N_d_2p_1;
    inputFile2 >> N_alpha_2s_2p_2 >> N_d_2s_2 >> N_d_2p_2;
    inputFile2 >> N_alpha_2s_2p_3 >> N_d_2s_3 >> N_d_2p_3;
    inputFile2.close();

    ifstream inputFile3(O_STO3G);
    inputFile3 >> O_alpha_2s_2p_1 >> O_d_2s_1 >> O_d_2p_1;
    inputFile3 >> O_alpha_2s_2p_2 >> O_d_2s_2 >> O_d_2p_2;
    inputFile3 >> O_alpha_2s_2p_3 >> O_d_2s_3 >> O_d_2p_3;
    inputFile3.close();

    ifstream inputFile4(F_STO3G);
    inputFile4 >> F_alpha_2s_2p_1 >> F_d_2s_1 >> F_d_2p_1;
    inputFile4 >> F_alpha_2s_2p_2 >> F_d_2s_2 >> F_d_2p_2;
    inputFile4 >> F_alpha_2s_2p_3 >> F_d_2s_3 >> F_d_2p_3;
    inputFile4.close();
    vector<vec> data;
    data.push_back(vec{H_alpha_1s_1,H_alpha_1s_2,H_alpha_1s_3});
    data.push_back(vec{H_d_1s_1,H_d_1s_2,H_d_1s_3});

    data.push_back(vec{C_alpha_2s_2p_1,C_alpha_2s_2p_2,C_alpha_2s_2p_3});
    data.push_back(vec{C_d_2s_1       ,C_d_2s_2       ,C_d_2s_3       });
    data.push_back(vec{C_d_2p_1       ,C_d_2p_2       ,C_d_2p_3       });

    data.push_back(vec{N_alpha_2s_2p_1,N_alpha_2s_2p_2,N_alpha_2s_2p_3});
    data.push_back(vec{N_d_2s_1       ,N_d_2s_2       ,N_d_2s_3       });
    data.push_back(vec{N_d_2p_1       ,N_d_2p_2       ,N_d_2p_3       });

    data.push_back(vec{O_alpha_2s_2p_1,O_alpha_2s_2p_2,O_alpha_2s_2p_3});
    data.push_back(vec{O_d_2s_1       ,O_d_2s_2       ,O_d_2s_3       });
    data.push_back(vec{O_d_2p_1       ,O_d_2p_2       ,O_d_2p_3       });

    data.push_back(vec{F_alpha_2s_2p_1,F_alpha_2s_2p_2,F_alpha_2s_2p_3});
    data.push_back(vec{F_d_2s_1       ,F_d_2s_2       ,F_d_2s_3       });
    data.push_back(vec{F_d_2p_1       ,F_d_2p_2       ,F_d_2p_3       });
    return data;
}

// Select exp, contraction, quantum num, and coord constants for a given basis function
tuple<vec,vec,vec,vec,vec> FindConsts(vector<vec> basis_data, int atom, vec R_center,string orbital){
    // Exponent Data
    vec H_alpha_1s    = basis_data[0];
    vec X_alpha_2s_2p;
    // Contraction Coefficient Data
    vec H_d_1s = basis_data[1];
    vec X_d_2s;
    vec X_d_2p;
    if (atom == 6){
        X_alpha_2s_2p = basis_data[2];
        X_d_2s = basis_data[3];
        X_d_2p = basis_data[4];
    } else if (atom == 7){
        X_alpha_2s_2p = basis_data[5];
        X_d_2s = basis_data[6];
        X_d_2p = basis_data[7];
    } else if (atom == 8){
        X_alpha_2s_2p = basis_data[8];
        X_d_2s = basis_data[9];
        X_d_2p = basis_data[10];
    } else if (atom == 9){
        X_alpha_2s_2p = basis_data[11];
        X_d_2s = basis_data[12];
        X_d_2p = basis_data[13];
    }

    // Quantum Numbers
    const vec lms_s = vec({0,0,0});
    const vec lms_px = vec({1,0,0});
    const vec lms_py = vec({0,1,0});
    const vec lms_pz = vec({0,0,1});

    // Initialize Output Vectors
    vec exponents;
    vec contraCoeffs;
    vec quantNums;
    if (atom == 1){
        exponents    = H_alpha_1s;
        contraCoeffs = H_d_1s;
        quantNums    = lms_s;
    } else if (atom == 6 || atom == 7 || atom == 8 || atom == 9){
        if        (orbital == "2s" ){
            exponents    = X_alpha_2s_2p ;
            contraCoeffs = X_d_2s;
            quantNums    = lms_s;
        } else if (orbital == "2px"){
            exponents    = X_alpha_2s_2p;
            contraCoeffs = X_d_2p;
            quantNums    = lms_px;
        } else if (orbital == "2py"){
            exponents    = X_alpha_2s_2p;
            contraCoeffs = X_d_2p;
            quantNums    = lms_py;
        } else if (orbital == "2pz"){
            exponents    = X_alpha_2s_2p;
            contraCoeffs = X_d_2p;
            quantNums    = lms_pz;
        }
    }
    // For each basis function, ωμ(r) (for μ = 1 · · · N ) there will be ...

    // (i) a center, R,
    double X = R_center(0);
    double Y = R_center(1);
    double Z = R_center(2);

    // (ii) 3 quantum numbers, (l , m, n),
    double l = quantNums(0);
    double m = quantNums(1);
    double n = quantNums(2);

    // and (iii) information about 3 primitive functions:
    // 3 exponents, αk,
    // 3 corresponding contraction coefficients, dk, and
    // 3 normalization constants, N_k_lmn
    int K = contraCoeffs.size();
    vec normConsts(K, fill::ones);
    for (int k = 0; k < K; k++){
        double alpha_k = exponents(k);
        double N_k_lmn = NormConst(X,Y,Z,alpha_k,l,m,n);
        normConsts(k) = N_k_lmn;
    }

    return make_tuple(R_center,quantNums,exponents,contraCoeffs,normConsts);
}

// Read in input file
tuple<int,vector<int>,int,int,vector<vector<double>>,vector<int>,vector<vector<double>>,vector<int>,vector<int>> read_input_file(string file_name){
    // Read in the coordinates, in the format: E X Y Z for each atom, where E is the element (handle at least H and C).

    // Read in input file
    ifstream inputFile(file_name);

    // Initialize vars
    int n_atoms;          // Initialize total number of atoms
    int charge;
    inputFile >> n_atoms >> charge; // Set total number of atoms and charge
    const double Bohr_A = 0.52917706; // angstroms in 1 bohr
    int a = 0; // Number carbons
    int b = 0;
    vector<vector<double>> xyz_list;       // Initialize list for atoms' xyz coordinates
    vector<int> atom_list;                 // Initialize list for atoms' identities
    vector<vector<double>> basis_xyz_list; // Initialize list for atoms' xyz coordinates
    vector<int> basis_atom_list;           // Initialize list for atoms' identities
    vector<int> basis_atom_nums;
    vector<int> Z;
    // Read in atom identity and xyz coordinates
    for (int i = 0; i < n_atoms; ++i) {  // Iterate through every atom
        int atom;
        double x, y, z;                    // Initialize atom identity and xyz coordinates
        inputFile >> atom >> x >> y >> z ; // Set atomic number/atom identity and xyz coordinates
        //x = x/Bohr_A;
        //y = y/Bohr_A;
        //z = z/Bohr_A;
        if (atom != 1 && atom != 6 && atom != 7 && atom != 8 && atom != 9) {      // If a given atom is not H or C, throw an error
            cerr << "Atom No." << i+1 << ": This atom is not a hydrogen, carbon, nitrogen, oxygen, fluorine!" << endl;
        } else if (atom == 6 || atom == 7 || atom == 8 || atom == 9){
            a = a + 1;
            basis_atom_list.insert(basis_atom_list.end(), {atom,atom,atom,atom});
            basis_xyz_list.insert(basis_xyz_list.end(), { {x, y, z},{x, y, z},{x, y, z},{x, y, z} });
            basis_atom_nums.insert(basis_atom_nums.end(),{i,i,i,i});
        } else if (atom == 1){
            b = b + 1;
            basis_atom_list.push_back(atom);
            basis_xyz_list.push_back({x, y, z});
            basis_atom_nums.push_back(i);
        }

        atom_list.push_back(atom);         // Append this atom's atomic number/atom identity to list
        xyz_list.push_back({x, y, z});     // Append this atom's xyz coordinates to list

        // ZA
        if (atom == 1){
            Z.push_back(1);
        } else if (atom == 6){
            Z.push_back(4);
        } else if (atom == 7){
            Z.push_back(5);
        } else if (atom == 8){
            Z.push_back(6);
        } else if (atom == 9){
            Z.push_back(7);
        }
    }
    inputFile.close();                     // Close the txt file

    tuple<int,vector<int>,int,int,vector<vector<double>>,vector<int>,vector<vector<double>>,vector<int>,vector<int>> output = make_tuple(n_atoms,Z,a,b,xyz_list,atom_list,basis_xyz_list,basis_atom_list,basis_atom_nums);
    return output;
}

tuple<vector<tuple<vec,vec,vec,vec,vec>>,vector<string>> build_basis(vector<vec> basis_data, int N, vector<int> basis_atom_list, vector<vector<double>> basis_xyz_list){
    //input basis_data, N, basis_atom_list, basis_xyz_list
    // Build a list of the basis functions, which are contracted gaussians,
    // Iterate through N AOs and construct basis functions
    vector<string> row2_orbital_bank = {"2s","2px","2py","2pz"};
    vector<tuple<vec,vec,vec,vec,vec>> basis_func_constants;
    vector<string> basis_orbital_list;


    for (int i = 0; i < N; i++){

        int atom = basis_atom_list[i]; // select atom
        vector<double> center = basis_xyz_list[i]; // find coordinates

        string orbital;
        if (atom == 1){
            orbital = "1s"; // If it's H, only 1s orbital
        } else if (atom == 6){ // If it's C, for each C, there's 2s, 2px, 2py, and 2pz
            orbital = row2_orbital_bank[0];
            row2_orbital_bank.erase(row2_orbital_bank.begin());
        }
        if (row2_orbital_bank.size() == 0){
            row2_orbital_bank = {"2s","2px","2py","2pz"};
        }
        basis_orbital_list.push_back(orbital); // Collect basis functions' orbitals + find their characteristic consts
        basis_func_constants.push_back(FindConsts(basis_data,atom, center, orbital));
    }
    return make_tuple(basis_func_constants,basis_orbital_list);
}

mat build_s(int N, vector<tuple<vec,vec,vec,vec,vec>> basis_func_constants){
    // Find all normalization constants for the basis functions
    mat All_Normalization_Constants(3,N,fill::zeros);

    // Run a loop over all your basis functions
    for (int i = 0; i < N; i++){
        // get the normalization constants for the 3 primitives that make up each basis function
        auto[R_center,quantNums,exponents,contraCoeffs,normConsts] = basis_func_constants[i];
        // and save them in an array.
        All_Normalization_Constants.col(i) = normConsts;
    }


    // Contracted overlap integral S
    mat S(N,N,fill::zeros);

    // Iterate over pairs of basis functions
    for (int mu = 0; mu < N; mu++){
        for (int nu = 0; nu < N; nu++){
            // For this given pair of basis functions
            // R_center,quantNums,exponents,contraCoeffs,normConsts for each atom in the pair
            auto[R_mu,lmn_mu,a_mu,d_mu,N_mu] = basis_func_constants[mu];
            auto[R_nu,lmn_nu,a_nu,d_nu,N_nu] = basis_func_constants[nu];

            // # of contracted gaussians
            int K = d_mu.size();
            int L = d_nu.size();

            double S_mu_nu = 0;

            for (int k = 0; k < K; k++){
                for (int l = 0; l < L; l++){
                    double d_k_mu     = d_mu(k);
                    double d_l_nu     = d_nu(l);
                    double N_k_mu     = N_mu(k);
                    double N_l_nu     = N_nu(l);
                    double alpha_k_mu = a_mu(k);
                    double alpha_l_nu = a_nu(l);

                    double Xm = R_mu(0);
                    double Ym = R_mu(1);
                    double Zm = R_mu(2);

                    double Xn = R_nu(0);
                    double Yn = R_nu(1);
                    double Zn = R_nu(2);

                    double lm = lmn_mu(0);
                    double mm = lmn_mu(1);
                    double nm = lmn_mu(2);

                    double ln = lmn_nu(0);
                    double mn = lmn_nu(1);
                    double nn = lmn_nu(2);

                    double S_k_l = 1; // Initialize SAB product
                    S_k_l *= SxABComp(Xm,Xn,alpha_k_mu,alpha_l_nu,lm,ln); // Compute Sxkl
                    S_k_l *= SxABComp(Ym,Yn,alpha_k_mu,alpha_l_nu,mm,mn); // Compute Sykl
                    S_k_l *= SxABComp(Zm,Zn,alpha_k_mu,alpha_l_nu,nm,nn); // Compute Szkl

                    S_mu_nu += d_k_mu * d_l_nu * N_k_mu * N_l_nu * S_k_l;
                }
            }
            S(mu,nu) = S_mu_nu;
        }
    }
    return S;
}

tuple<double,double,double,double,double,double,double,double>gamma_atom_constants(vector<vec> basis_data,vector<int> atom_list, vec RA, vec RB, int A, int B, int k, int kp, int l, int lp){
    vec d_A;
    vec d_B;
    vec alpha_A;
    vec alpha_B;

    int atom = atom_list[A];
    vec d;
    vec alpha;
    if (atom == 1){
        alpha = basis_data[0];
        d = basis_data[1];
    } else if (atom == 6){
        alpha = basis_data[2];
        d = basis_data[3];
    } else if (atom == 7){
        alpha = basis_data[5];
        d = basis_data[6];
    } else if (atom == 8){
        alpha = basis_data[8];
        d = basis_data[9];
    } else if (atom == 9){
        alpha = basis_data[11];
        d = basis_data[12];
    }
    d_A = d;
    alpha_A = alpha;

    atom = atom_list[B];
    if (atom == 1){
        alpha = basis_data[0];
        d = basis_data[1];
    } else if (atom == 6){
        alpha = basis_data[2];
        d = basis_data[3];
    } else if (atom == 7){
        alpha = basis_data[5];
        d = basis_data[6];
    } else if (atom == 8){
        alpha = basis_data[8];
        d = basis_data[9];
    } else if (atom == 9){
        alpha = basis_data[11];
        d = basis_data[12];
    }
    d_B = d;
    alpha_B = alpha;

    // Extract contraction coefficients
    double d_k_sA  = d_A(k);
    double d_kp_sA = d_A(kp);
    double d_l_sB  = d_B(l);
    double d_lp_sB = d_B(lp);
    cout << "k=" << k << ",d_k=" << d_k_sA << endl;
    cout << "kp=" << kp << ",d_kp=" << d_kp_sA << endl;
    cout << "l=" << l << ",d_l=" << d_l_sB << endl;
    cout << "lp=" << lp << ",d_lp=" << d_lp_sB << endl;
    // Extract exponents
    double alpha_k  = alpha_A(k);
    double alpha_kp = alpha_A(kp);
    double alpha_l  = alpha_B(l);
    double alpha_lp = alpha_B(lp);

    // Calculate Normalization Constants
    double N_k_sA  = NormConst(RA(0), RA(1), RA(2), alpha_k , 0, 0, 0);
    double N_kp_sA = NormConst(RA(0), RA(1), RA(2), alpha_kp, 0, 0, 0);
    double N_l_sB  = NormConst(RB(0), RB(1), RB(2), alpha_l , 0, 0, 0);
    double N_lp_sB = NormConst(RB(0), RB(1), RB(2), alpha_lp, 0, 0, 0);

    // Calculate d primes
    double dp_k_sA  = d_k_sA *N_k_sA;
    double dp_kp_sA = d_kp_sA*N_kp_sA;
    double dp_l_sB  = d_l_sB *N_l_sB;
    double dp_lp_sB = d_lp_sB*N_lp_sB;

    return make_tuple(alpha_k,alpha_kp,alpha_l,alpha_lp,dp_k_sA,dp_kp_sA,dp_l_sB,dp_lp_sB);
}

mat build_g(int n_atoms, vector<vector<double>> xyz_list,vector<vec> basis_data,vector<int> atom_list){
    // Contracted overlap integral S

    mat gamma(n_atoms,n_atoms,fill::zeros);

    for (int A = 0; A < n_atoms; A++){
        for (int B = 0; B < n_atoms; B++){
            //Define RA and RB
            vector<double> RA_raw = xyz_list[A];
            vector<double> RB_raw = xyz_list[B];
            vec RA = vec{RA_raw[0],RA_raw[1],RA_raw[2]};
            vec RB = vec{RB_raw[0],RB_raw[1],RB_raw[2]};

            double g = 0;

            for             (int k = 0; k < 3; k++){
                for         (int kp = 0; kp < 3; kp++){
                    for     (int l = 0; l < 3; l++){
                        for (int lp = 0; lp < 3; lp++){
                            double summand;

                            auto[alpha_k,alpha_kp,alpha_l,alpha_lp,dp_k_sA,dp_kp_sA,dp_l_sB,dp_lp_sB]= gamma_atom_constants(basis_data,atom_list,RA,RB,A,B,k,kp,l,lp);

                            double sigma_A = pow(alpha_k+alpha_kp,-1);
                            double sigma_B = pow(alpha_l+alpha_lp,-1);

                            double UA = ((M_PI*sigma_A),3/2);
                            double UB = ((M_PI*sigma_B),3/2);

                            double V2 = pow(sigma_A+sigma_B,-1);

                            double T = V2*dot(RA-RB,RA-RB);
                            double oo;
                            if (A==B){
                                oo = UA*UB*sqrt(2*V2)*sqrt(2/M_PI);
                            } else {
                                oo = UA*UB*sqrt(1/(dot(RA-RB,RA-RB)))*erf(sqrt(T));
                            }
                            summand = dp_k_sA*dp_kp_sA*dp_l_sB*dp_lp_sB*oo;
                            g += summand;
                        }
                    }
                }
            }
            gamma(A,B) = g;
        }
    }
    return gamma;
}


int main(int argc, char* argv[]) {

    // Program inputs
    //if (argc !=2)
    //{
    //    printf("Usage: hw3 <filename>, for example hw3 example.txt\n");
    //    return EXIT_FAILURE;
    //}
    //string file_name = argv[1];
    string file_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW4/sample_input/H2.txt";

    string H_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW4/basis/H_STO3G.txt";
    string C_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW4/basis/C_STO3G.txt";
    string N_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW4/basis/N_STO3G.txt";
    string O_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW4/basis/O_STO3G.txt";
    string F_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW4/basis/F_STO3G.txt";

    // Question 1
    // Read input file
    auto[n_atoms,Z,a,b,xyz_list,atom_list,basis_xyz_list,basis_atom_list,basis_atom_nums] = read_input_file(file_name);

    // Evaluate the number of basis functions, N from the molecular formula, Ca Hb , where
    // the relation is N = 4a + b. Your matrices, such as S, H, X, etc, will be N × N , so this will
    // enable you to define them.
    int N = 4*a+b;

    // Evaluate the number of electrons 2n = 4a + b. Throw an error if the number of electron
    // pairs n = 2a +b/2 is not an integer. Knowing n is necessary to evaluate the energy later.
    int n;
    if (!is_integer(2*a+b/2)){
        cerr << "The number of electron pairs n = 2a +b/2 is not an integer!" << endl;
    } else {
        n = 2*a + b/2;
    }
    // Build Basis Functions
    vector<vec> basis_data = ProcessBasisData(H_path_name,C_path_name, N_path_name, O_path_name, F_path_name);
    auto[basis_func_constants,basis_orbital_list] = build_basis(basis_data,N,basis_atom_list,basis_xyz_list);


    mat IAb = mat("7.176 14.051 19.316 25.390 32.272; NaN 5.572 7.275 9.111 11.080; 9 21 25 31 39");

    mat gamma = build_g(n_atoms,xyz_list,basis_data,atom_list);
    gamma.print("gamma");

    // Build Contracted Overlap Integral Matrix
    mat S = build_s(N,basis_func_constants);
    S.print("Overlap, S: Contracted overlap integral, overlap matrix:");

}
