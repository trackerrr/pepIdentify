//
//  main.cpp
//  CS482 A4
//
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

// aa table
float aa_mass(char aa) {
    float mass = 0;
    if (aa == 'A') {
        mass = 71.03711;
    }
    else if (aa == 'R') {
        mass = 156.10111;
    }
    else if (aa == 'N') {
        mass = 114.04293;
    }
    else if (aa == 'D') {
        mass = 115.02694;
    }
    else if (aa == 'C') {
        mass = 160.03065;
    }
    else if (aa == 'E') {
        mass = 129.04259;
    }
    else if (aa == 'Q') {
        mass = 128.05858;
    }
    else if (aa == 'G') {
        mass = 57.02146;
    }
    else if (aa == 'H') {
        mass = 137.05891;
    }
    else if (aa == 'I') {
        mass = 113.08406;
    }
    else if (aa == 'L') {
        mass = 113.08406;
    }
    else if (aa == 'K') {
        mass = 128.09496;
    }
    else if (aa == 'M') {
        mass = 131.04049;
    }
    else if (aa == 'F') {
        mass = 147.06841;
    }
    else if (aa == 'P') {
        mass = 97.05276;
    }
    else if (aa == 'S') {
        mass = 87.03203;
    }
    else if (aa == 'T') {
        mass = 101.04768;
    }
    else if (aa == 'W') {
        mass = 186.07931;
    }
    else if (aa == 'Y') {
        mass = 163.06333;
    }
    else if (aa == 'V') {
        mass = 99.06841;
    }
    return mass;
}
float pep_mass(string pep){
    float mass = 0;
    for (int i = 0; i < pep.size(); i++) {
        mass += aa_mass(pep[i]);
    }
    return mass+18.01;
}
// may not be used
struct pep_mass_rec {
    int combined_num, start_pos;
    float pep_mass;
    pep_mass_rec(int combined_num, int start_pos, float pep_mass){
        this->combined_num=combined_num;
        this->start_pos=start_pos;
        this->pep_mass=pep_mass;
    }
};
struct spectrum {
    // format:
    // id: (index) 0, 1...
    // peptide: TGAGNPVGDKLNVITVGPRGPLLV...
    // protein: (first 10 characters) >Q15843ups
    // m_vs_z: (PEPMASS=) 400.6561
    // z: (CHARGE=) 3+
    // num: scans number
    int id, z;
    string peptide, protein;
    float m_vs_z, spec_mass, score;
    float tallest_peak;
    vector<pair<float, float> > peak;
    
    spectrum(int id, int z, float m_vs_z){
        this->id = id;
        this->z = z;
        this->m_vs_z = m_vs_z;
        this->spec_mass = m_vs_z * z - 1.0073 * z;
        this->tallest_peak = 0;
        //cout<<mass<<endl;
    }
    
};
struct protein {
    string name, seq;
    vector<pair<int, int> > digested_pt; // digested point, starting and length
    int pep_num;
    float *pep_mass_records;
    // y-ions
    // pep1: ASDFG y-ions
    // pep2: EUWOHVB y-ions
    vector<float> *y_ion_mass;
    
    protein(string name, string seq){
        this->name = name;
        this->seq = seq;
        
    }
    void digest(){
        int start = 0; // start position
        for (int i = 0; i < seq.size(); i++) {
            if ((seq[i]=='R' || seq[i]=='K') && (i!=seq.size()-1 && seq[i+1]!='P')) {
                digested_pt.push_back(make_pair(start, i-start+1));
                start = i+1;
            }
            else if (i == seq.size()-1){
                digested_pt.push_back(make_pair(start, i-start+1));
            }
        }
        pep_num = digested_pt.size();
        pep_mass_records = new float [pep_num];
        y_ion_mass = new vector<float> [pep_num];
        // records each possible mass
        // total: n peptide pieces
        // first n original peptides mass
        for (int i = 0; i < pep_num; i++) {
            // get the position of digested_pep[i]
            int s = digested_pt.at(i).first;
            int l = digested_pt.at(i).second;
            string cur_pep = seq.substr(s,l);
            float cur_pep_mass = pep_mass(cur_pep);
            pep_mass_records[i] = cur_pep_mass;
            //cout<<s<<" "<<l<<endl;
            //cout<<cur_pep<<" "<<cur_pep_mass<<endl;
        }
    }
    void find(vector<spectrum> *spec_list) {
        // find each peptide in the protein
        // j is the order of peptides
        for (int j = 0; j < pep_num; j++) {
            float pep_m = pep_mass_records[j] - 18.01 + 19.0178; // without water, with y-ion add on
            int s = digested_pt.at(j).first; // starting in seq
            int l = digested_pt.at(j).second;
            
            // search in all spectra
            // k is the spectrum position
            for (int k = 0; k < spec_list->size(); k++) {
                spectrum spec = spec_list->at(k);
                float error = spec.spec_mass - pep_mass_records[j];
                
                // we find a match
                if (abs(error) <= 0.1 ) {
                    //cout<<"pep id: "<<j<< ", spec id: " <<spec.id <<endl;
                    //cout << spec.spec_mass <<endl;
                    
                    // y-ion mass[j]
                    // if we have NOT calculate this peptides' y-ions
                    if (y_ion_mass[j].size()==0) {
                        
                        //cout<<"pep_num: "<<j <<endl;
                        //cout<<"spec: "<<k<<endl;
                        //cout<<"pep_m: "<<pep_m<<endl;
                        //cout << seq.substr(s,l)<<endl;
                        float cur_mass = pep_m;
                        for (int pos = 0; pos < l-1; pos++) {
                            // pos is the one to be deleted
                            cur_mass -= aa_mass(seq[s+pos]);
                            y_ion_mass[j].push_back(cur_mass);
                            //cout << "y ion mass: "<<cur_mass<<endl;
                        }
                    }
                    
                    // calculate score for a spectrum
                    float score = 0;
                    for (int x = 0; x < y_ion_mass[j].size(); x++) {
                        // search for each y-ion
                        float y_m = y_ion_mass[j].at(x);
                        float tall_peak_intensity = 0;
                        // go thru spectrum peak list
                        for (int w = 0; w < spec.peak.size();w++) {
                            float peak = spec.peak.at(w).first;
                            float intensity = spec.peak.at(w).second;
                            if (abs(peak - y_m)<=0.5) {
                                if (intensity > tall_peak_intensity) {
                                    tall_peak_intensity = intensity;
                                    //cout<<peak << " " << y_m<<endl;
                                }
                            }
                            if (peak - y_m>0.5){
                                break;
                            }
                        }
                        if(log10(100*tall_peak_intensity/spec.tallest_peak)>0){
                            score += log10(100*tall_peak_intensity/spec.tallest_peak);
                            //cout<<score<<endl;
                        }
                        
                    }
                    if (score > spec.score) {
                        // update in spec list
                        spec_list->at(k).score = score;
                        spec_list->at(k).protein = name;
                        spec_list->at(k).peptide = seq.substr(s,l);
                        //cout << "spec: " << k <<endl;
                        //cout << "protein " << spec_list->at(k).protein <<endl;
                        //cout << "peptide " << spec_list->at(k).peptide <<endl;
                        //cout <<  "score: " << spec_list->at(k).score <<endl;
                    }
                }
            }
            
        }
        
    }
};



int main(int argc, const char * argv[]) {
    ifstream mgf,fasta;
    mgf.open(argv[1]);
    fasta.open(argv[2]);
    vector<spectrum> *spec_list = new vector<spectrum>;
    vector<protein> *pro_list = new vector<protein>;
    
    // read in mgf
    string line;
    int id = 0;
    int z;
    float m_vs_z;
    string::size_type sz; // alias of m_vs_z
    
    while(!mgf.eof()) {
        getline(mgf,line);
        if(line != "BEGIN IONS"){
            break;
        }
        // discard "BEGIN IONS"
        // discard index
        getline(mgf,line);
        // read in m/z
        getline(mgf,line);
        line = line.substr(8);
        m_vs_z = stod(line,&sz);
        
        // read in charge
        getline(mgf,line);
        char c = line[7];
        z = c - 48;
        // skip scan
        getline(mgf,line);
        getline(mgf,line);
        
        spectrum spec(id, z, m_vs_z);
        getline(mgf,line);
        while (line != "END IONS") {
            float fst, sec;
            fst = stod(line,&sz);
            sec = stod(line.substr(sz));
            spec.peak.push_back(make_pair(fst, sec));
            // find tall peak
            if (sec>spec.tallest_peak) {
                spec.tallest_peak=sec;
            }
            getline(mgf,line);
        }
        //cout <<spec.tallest_peak<<endl;break;
        spec_list->push_back(spec);
        // discard new line
        getline(mgf,line);
        id++;
    }
    mgf.close();
    
    // read in fasta
    string name, seq;
    while (!fasta.eof()) {
        getline(fasta, line);
        name = line.substr(0,10);
        getline(fasta, seq);
        struct protein pro(name,seq);
        pro_list->push_back(pro);
    }
    fasta.close();
    
    // find for all proteins
    for (int o = 0; o < pro_list->size(); o++) {
        pro_list->at(o).digest();
        pro_list->at(o).find(spec_list);
    }
    
    // output
    for (int i = 0; i < spec_list->size(); i++) {
        spectrum spec = spec_list->at(i);
        if (spec.peptide != "") {
            cout << spec.id << " " << spec.m_vs_z << " " << spec.z << " " << spec.peptide << " " << spec.protein << " " << spec.score <<endl;
        }
        
    }
    
}























