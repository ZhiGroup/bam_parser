//
//  HL.h
//  SL_test
//
//  Created by Degui Zhi on 7/25/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef SL_test_HL_h
#define SL_test_HL_h

class CompCol_Mat_double;
class CompRow_Mat_double;
class Coord_Mat_double;

class HL {
    
private:
    Coord_Mat_double      seq_;				// the actual base-calls
    CompCol_Mat_double    seq_Col_;			// the actual base-calls
    CompRow_Mat_double    seq_Row_;			// the actual base-calls
    Coord_Mat_double      ql_;				// the quality values
    CompCol_Mat_double    ql_Col_;			// the quality values
    CompRow_Mat_double    ql_Row_;			// the quality values
	
    VECTOR_int		read_;			// read index
    VECTOR_int		pos_;			// row_ind (nz_ elements)
    char*			chr_;    
    
    int nz_;                   // number of nonzeros, to be made equal for the data matrices
    int dim_[2];               // number of rows, cols, to be made equal for the data matrices
    
public:
    HL(void);
    HL(const Coord_Mat_double &S, const Coord_Mat_double &Q);
	HL(const char *bamfile, const char *markerfile); // construct HL from BAM file over predefined markers
    ~Coord_Mat_double() {};
    
    /*******************************/
    /*  Access and info functions  */
    /*******************************/
    
    double&      val(int i) { return val_(i); }
    int&         row_ind(int i) { return rowind_(i); }
    int&         col_ind(int i) { return colind_(i);}
    
    const double&      val(int i) const { return val_(i); }
    const int&         row_ind(int i) const { return rowind_(i); }
    const int&         col_ind(int i) const { return colind_(i);}
    
    int          dim(int i) const {return dim_[i];};
    int          size(int i) const {return dim_[i];};
    int          NumNonzeros() const {return nz_;};
	
	int		single_read(int i); // return the indices of single-site reads
	int		double_read(int i); // return the indices of double-site reads
	
    /*******************************/
    /* I/O  */
    /*******************************/
    
	void load(const char *filename);
	void write(const char *filename);	
	void write_VCF(const char *vcffilename);	
	
    /*******************************/
    /* GL and HL calculation  */
    /*******************************/
	
	VECTOR_double calculate_GL (); // for all, of size 3 * num_sites
	VECTOR_double calculate_GL (VECTOR_int pos_ind); // for a subset
	VECTOR_double calculate_GL (VECTOR_int read_ind); // for a subset
	VECTOR_double calculate_GL (VECTOR_int pos_ind, VECTOR_int read_ind); // for a subset
	
	VECTOR_double calculate_HL (char oddness); // for all, of size 10 * num_sites
	VECTOR_double calculate_HL (VECTOR_int pos_ind, char oddness); // for a subset	
	VECTOR_double calculate_HL (VECTOR_int read_ind, char oddness); // for a subset	
	VECTOR_double calculate_HL (VECTOR_int pos_ind, VECTOR_int read_ind, char oddness); // for a subset	
	
    /***********************************/
    /*  Operations         */
    /***********************************/
	HL& subset(const HL &H, VECTOR_int pos_ind, VECTOR_int read_ind);
	HL& partition();	// into GL, HL0, and HL1
    	
    /***********************************/
    /*  General access function (slow), may not needed */ 
    /***********************************/
    
    double       operator() (int i, int j) const;        
    double&      set(int i, int j);
        
};





#endif
