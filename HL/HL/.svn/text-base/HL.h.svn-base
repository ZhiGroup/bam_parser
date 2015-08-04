// To Dos:
/*
 * do GL
 * do HL
 * check paired-end reads
 * memory allocation
 * skip indels
 * do partition
 * put ref_id in hash 
 */

//
//  HL.h
//  HL
//
//  Created by Degui Zhi on 7/26/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef HL_HL_h
#define HL_HL_h

#include "vecdefs.h"
#include VECTOR_H
#include "compcol_double.h"
#include "comprow_double.h"
#include "coord_double.h"
#include <string.h>
#include <stdlib.h>
#include <vector.h>
#include "MyHash.h"

using namespace std;

class CompCol_Mat_double;
class CompRow_Mat_double;
class Coord_Mat_double;

class HL {
    
private:
	
    vector<int>		read_;			// read index, // useless so far
    // read_ are supposedly sorted by starting positions
    vector<string>  read_name_;     // read name from BAM, need to check if read pairs are handled well
    vector<int>            pos_;			// row_ind
    char*			chr_;    
    vector<string>  ref_;
    vector<string>  alt_;
    
    int nz_;                   // number of nonzeros, to be made equal for the data matrices
    int dim_[2];               // number of rows, cols, to be made equal for the data matrices
    
public:
    HL(void);
    HL(const Coord_Mat_double &S, const Coord_Mat_double &Q);
    ~HL();
    
    // current practice: only save the Coord_Mat_double, and generate the rest as needed
    Coord_Mat_double      seq_;				// the actual base-calls
    CompCol_Mat_double    seq_Col_;			// the actual base-calls
    CompRow_Mat_double    seq_Row_;			// the actual base-calls
    Coord_Mat_double      s01_;				// the sequence in 01, ignored non ref/alt
    CompCol_Mat_double    s01_Col_;			// the sequence in 01
    CompRow_Mat_double    s01_Row_;			// the sequence in 01
    Coord_Mat_double      ql_;				// the quality values
    CompCol_Mat_double    ql_Col_;			// the quality values
    CompRow_Mat_double    ql_Row_;			// the quality values
    Coord_Mat_double      q01_;				// the quality values, ignored non ref/alt
    CompCol_Mat_double    q01_Col_;			// the quality values
    CompRow_Mat_double    q01_Row_;			// the quality values
    
    // more efficiently would be to record the values seq_.val() and 
    
	int parse_bam(const char *bamfile, const char *markerfile); // construct HL from BAM file over predefined markers
	int parse_bam_vcf(const char *bamfile, const char *vcfmarkerfile); // construct HL from BAM file over predefined markers
	int read_bam_vcf(const char *bamfile, const char *vcfmarkerfile, const char *hlfile, const char *readfile, const char *sitefile);
    int seq_to_s01();
    
    /*******************************/
    /*  Access and info functions  */
    /*******************************/
        
    int          dim(int i) const {return dim_[i];};
    int          size(int i) const {return dim_[i];};
    int          NumNonzeros() const {return nz_;};
	
	int		get_single_reads(vector<int>& read_subset); // return the indices of single-site reads
	int		get_double_reads(vector<int>& read_subset); // return the indices of double-site reads
	int		get_long_reads(vector<int>& read_subset); // return the indices of reads spanning >2 sites
	vector<int>	merge_read_sets(const vector<int>& read_subset1, const vector<int>& read_subset2); // merge two sorted indices
	
    /*******************************/
    /* I/O  */
    /*******************************/
    
	void load(const char *filename); 
	void write(const char *filename, const char *style);	//done
	void write_VCF(const char *vcffilename);	
    void write_sites(FILE *out_file, const char *style);    // done
	void write_reads(FILE *out_file);                       // done
	void write_counts(const int person_id, const char *outfilename);                       // done
	void write_jumps(const int person_id, const char *of1, const char *of2);                       // done
	void write_matrix(FILE *out_file, const char *style);   // done
    void write_hl_by_read(const int person_id, const char *outfilename, const int use_quality);
	
    /*******************************/
    /* GL and HL calculation  */
    /*******************************/
	
	double** calculate_GL (); // for all, of size 3 * num_sites
	VECTOR_double calculate_GL_pos (vector<int> pos_ind); // for a subset
	VECTOR_double calculate_GL_read (vector<int> read_ind); // for a subset
	VECTOR_double calculate_GL (vector<int> pos_ind, vector<int> read_ind); // for a subset
	
	VECTOR_double calculate_HL (char oddness); // for all, of size 10 * num_sites
	VECTOR_double calculate_HL_pos (vector<int> pos_ind, char oddness); // for a subset	
	VECTOR_double calculate_HL_read (vector<int> read_ind, char oddness); // for a subset	
	VECTOR_double calculate_HL (vector<int> pos_ind, vector<int> read_ind, char oddness); // for a subset	
	
    /***********************************/
    /*  Operations         */
    /***********************************/
	HL subset(vector<int> pos_ind, vector<int> read_ind);
	HL subset_by_sites(vector<int> pos_ind);
	HL subset_by_reads(int * read_ind, int k);          // done
    HL subset_by_reads(vector<int>& read_ind, int M);
	int partition(HL& even, HL& odd);	// into even-split and odd-split
	int split(HL& h, int oddness);	// into even-split and odd-split, return number of reads
    
    /***********************************/
    /*  General access function (slow), may not needed */ 
    /***********************************/
    
    double       operator() (int i, int j) const;        
    double&      set(int i, int j);
    
};

#endif
