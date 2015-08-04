//
//  HL.cpp
//  HL
//
//  Created by Degui Zhi on 8/1/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "HL.h"
#include "sam.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream.h>
#include <string.h>
#include <hash_map.h>
#include <vector.h>
#include <algorithm>
#include <tr1/unordered_map>

using namespace std;

#include "compcol_double.h"
#include "comprow_double.h"
#include "coord_double.h"
#include "ilupre_double.h"
#include "icpre_double.h"
#include "iotext_double.h"


HL::HL()
{
    // read_ = NULL;			// read index
    // pos_=NULL;			// row_ind (nz_ elements)
    // chr_=NULL;        
    nz_=0;                   // number of nonzeros, to be made equal for the data matrices
    dim_[0] = 0; dim_[1] = 0;           // number of rows, cols, to be made equal for the data matrices
    
    /*
     seq_ = NULL;				// the actual base-calls
     seq_Col_ = NULL;			// the actual base-calls
     seq_Row_ = NULL;			// the actual base-calls
     ql_ = NULL;				// the quality values
     ql_Col_ = NULL;			// the quality values
     ql_Row_ = NULL;			// the quality values  
     */
}

HL::~HL()
{/*
  delete [] read_;
  delete [] pos_;
  delete [] chr_;
  delete [] seq_;
  delete [] seq_Col_;
  delete [] seq_Row_;
  delete [] ql_;
  delete [] ql_Row_;
  delete [] ql_Col_;
  */
}

struct eqstr
{
    bool operator()(const char* s1, const char* s2) const
    {
        return strcmp(s1, s2) == 0;
    }
};

// std::tr1::unordered_map<char*, int, eqstr> read_ind;
std::tr1::unordered_map<string, int> read_ind;
std::tr1::unordered_map<int, int> pos_ind;
// hash_map<const char*, int, hash<const char*>, eqstr > read_ind;
// hash_map<const int, int, hash<const char*>, eqstr > pos_ind;

int read_counter=0;
int pos_ptr = 0;
int site_id = -1;
int maxcol = 0;

// need memory allocation here
int ptr = 0;
vector<int> site_i;
vector<int> read_j;
vector<double> nuc;
vector<double> qual;
vector<string> read_name;
vector<int> pos_name;

int old_pos = -1;
int pos_ind1 = 0;


typedef struct {
    int beg, end;
    samfile_t *in;
    HL *hl;
    FILE *hlfile, *readfile, *sitefile;
} tmpstruct_t;

/*
 void add_memory (int ptr)
{
    site_i.reserve(ptr);
    read_j.reserve(ptr);
    nuc.reserve(ptr);
    qual.reserve(ptr);
    read_name.reserve(ptr);
    refs.reserve(ptr);
    alts.reserve(ptr);
    pos_name.reserve(ptr);
}
*/

// parse one line of VCF file
void vcf_parse_line(bam_header_t *header,  char *p, int *ref_id, int *pos, char *ref, char *alt)
{
    char *pch;
    pch = strtok (p,"\t\n");
    int item_id = 0;
//    strcpy(pch,p);
    while (pch != NULL)
    {
        item_id ++;
        if (item_id > 5) break;
        if (item_id > 5) break;
        if (item_id == 1) {
            *ref_id = -1;
            for (int i=0; i<header->n_targets; i++)
            {
                if (strcmp(pch, header->target_name[i])==0)
                {
                    *ref_id = i;
                    break;
                }
            }
            if (*ref_id==-1)
                cerr << "mismatched chr names in VCF and BAM" << endl;
        }
        // *ref_id = 0;
        if (item_id == 2) *pos = atoi(pch);
        // if (item_id == 3) continue;
        if (item_id == 4) strcpy(ref, pch);
        if (item_id == 5) strcpy(alt, pch);
        pch = strtok (NULL, "\t\n");
    }
}

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
    bam_plbuf_t *buf = (bam_plbuf_t*)data;
    bam_plbuf_push(b, buf);
    return 0;
}

static int pileup_func_novcf(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    int c,i,q;
    int old_ptr=0;
    char *chr;
    site_id ++;
    tmpstruct_t *tmp = (tmpstruct_t*)data;
    if ((int)pos >= tmp->beg && (int)pos < tmp->end) {
        chr  = tmp->in->header->target_name[tid];
        // printf("%s\t%d\t%d\t%d", chr, pos + 1, n, old_pos);
        //        if ((int)pos > maxcol) maxcol = pos;
        
        for (i = 0; i < n; ++i) {
            const bam_pileup1_t *p = pl + i;
            c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
            char *name = bam1_qname(p->b);
            q = bam1_qual(p->b)[p->qpos];
            
            // put pos,name,c,q into HL data structure
            
            if (read_ind.count(name)>0)
            {}
            else {
                // printf("new read %d:%s(%c,%d) \n", read_ind[name], name, c, q);
                
                read_name.push_back( name );
                read_ind[name] = ++read_counter;
                //                read_ind.insert( std::make_pair(name, ++read_counter) );
            }
            // printf("read_ind[%s]=%d @%d=(%c,%d) ", name, read_ind[name], pos, c, q);
            
            //if ((int)pos > old_pos) {
            old_pos = pos;
            pos_name.push_back( pos );
            pos_ind1 ++;
            maxcol=pos_ind1;
            
            // }
            site_i.push_back( pos_ind1-1);
            read_j.push_back( read_ind[name]-1);
            nuc.push_back(  c);
            qual.push_back( q );
            // cout << "size of qual = " << qual.size() << endl;
            ptr++;
        }
        // printf("\n");
        
        // add memory
        if ((ptr % 100000) < (old_ptr % 100000))
        {
            cout << ptr << "non zero entries loaded" << endl;
            // add_memory(ptr);
        }
        old_ptr = ptr;
    }
    return 0;
}


// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    int c,i,q;
    int old_ptr=0;
    char *chr;
    site_id ++;
    tmpstruct_t *tmp = (tmpstruct_t*)data;
    if ((int)pos >= tmp->beg && (int)pos < tmp->end) {
        chr  = tmp->in->header->target_name[tid];
        // printf("%s\t%d\t%d\t%d", chr, pos + 1, n, old_pos);
        //        if ((int)pos > maxcol) maxcol = pos;
        
        for (i = 0; i < n; ++i) {
            const bam_pileup1_t *p = pl + i;
            c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
            char *name = bam1_qname(p->b);
            q = bam1_qual(p->b)[p->qpos];
            
            // put pos,name,c,q into HL data structure
            
            if (read_ind.count(name)>0)
            {}
            else {
                // printf("new read %d:%s(%c,%d) \n", read_ind[name], name, c, q);
                
                read_name.push_back( name );
                read_ind[name] = ++read_counter;
                //                read_ind.insert( std::make_pair(name, ++read_counter) );
            }
            // printf("read_ind[%s]=%d @%d=(%c,%d) ", name, read_ind[name], pos, c, q);
            
            //if ((int)pos > old_pos) {
            // old_pos = pos;
            // pos_name.push_back( pos );
            // pos_ind ++;
            // maxcol=pos_ind;
            
            // }
            site_i.push_back( pos_ptr);
            read_j.push_back( read_ind[name]-1);
            nuc.push_back(  c);
            qual.push_back( q );
            // cout << "size of qual = " << qual.size() << endl;
            ptr++;
        }
        // printf("\n");
        
        // add memory
        if ((ptr % 100000) < (old_ptr % 100000))
        {
            cout << ptr << "non zero entries loaded" << endl;
            // add_memory(ptr);
        }
        old_ptr = ptr;
    }
    return 0;
}


static int pileup_func_nosave(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i,q;
    char c;
    int old_ptr=0;
    char *chr;
    site_id ++;
    tmpstruct_t *tmp = (tmpstruct_t*)data;
    if ((int)pos >= tmp->beg && (int)pos < tmp->end) {
        chr  = tmp->in->header->target_name[tid];
        for (i = 0; i < n; ++i) {
            const bam_pileup1_t *p = pl + i;
            c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
            char *name = bam1_qname(p->b);
            q = bam1_qual(p->b)[p->qpos];
            
            // put pos,name,c,q into HL data structure
            
            if (read_ind.count(name)>0)
            {}
            else {
                // printf("new read %d:%s(%c,%d) \n", read_ind[name], name, c, q);
                
                read_name.push_back( name );
                read_ind[name] = ++read_counter;
                //                read_ind.insert( std::make_pair(name, ++read_counter) );
            }
            // printf("read_ind[%s]=%d @%d=(%c,%d) ", name, read_ind[name], pos, c, q);
            
            if ((int)pos > old_pos) {
            old_pos = pos;
            pos_name.push_back( pos );
            pos_ind1 ++;
            maxcol=pos_ind1;
            
            }
            cout << pos_ind1 << "\t";
            cout << read_ind[name] << "\t";
            cout << c << "\t";
            cout << q << "\n";
  
            site_i.push_back( pos_ind1-1);
            read_j.push_back( read_ind[name]-1);
            nuc.push_back(  c);
            qual.push_back( q );
            // cout << "size of qual = " << qual.size() << endl;
            ptr++;
        }
        // printf("\n");
        
        old_ptr = ptr;
    }
    return 0;
}


int HL::parse_bam_vcf(const char *bamfile, const char *vcfmarkerfile)
// construct HL from BAM file 
// over predefined markers as defined by VCF-formatted file
// VCF file is sorted
// only VCF chr pos ref alt are used
{
	tmpstruct_t tmp;
	tmp.beg = 0; tmp.end = 0x7fffffff;
	tmp.in = samopen(bamfile, "rb", 0);
    vector<string> refs;
    vector<string> alts;
    vector<int> poss;
   
    if (tmp.in == 0) {
		fprintf(stderr, "vcf Fail to open BAM file %s\n", bamfile);
		return 1;
	}
    int ref;
    int beg;
    char *ref_allele = new char[1000];
    char *alt_allele = new char[1000];
    bam_index_t *idx;
    bam_plbuf_t *buf;
    idx = bam_index_load(bamfile); // load BAM index
    if (idx == 0) {
        fprintf(stderr, "BAM indexing file is not available.\n");
        return 1;
    }
    
	char *p = new char [100000]; // sometimes may not be big enough
					// Kui Zhang Modified
    ifstream myfile;
    myfile.open(vcfmarkerfile);
    string line;
    if (myfile.is_open())
    {
        while ( myfile.good() )
        {
            getline (myfile,line);
            if (line[0] == '#')
                continue;    
            if ( line.size()==0 )
                break;
            // cout << "read line: " << line << endl;
            // char *p=new char[line.size()]; // Kui Zhang Modified 
            strcpy(p,line.c_str());
            vcf_parse_line(tmp.in->header, p, &ref, &beg, ref_allele, alt_allele); // parse a variant
            tmp.beg = beg-1;
            tmp.end = beg;
            // cout << "ref=" << ref << "; beg=" << tmp.beg << endl;
            if (ref < 0) {
                fprintf(stderr, "Invalid region in file %s\n", vcfmarkerfile);
                return 1;
            }
            // ugly global pos_ind
//            refs[pos_ind] = ref_allele;
//            alts[pos_ind] = alt_allele;
            poss.push_back(beg);
            refs.push_back(ref_allele);
            alts.push_back(alt_allele);
            pos_ind[beg] = pos_ptr;
            //pos_ptr ++;
            // cout << "pos_ptr=" << pos_ptr << "; pos_ind.size()=" << pos_ind.size() << endl;
            
            buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
            bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
            bam_plbuf_push(0, buf); // finalize pileup
            pos_ptr ++;
        }
        myfile.close();
    }
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);
	
	samclose(tmp.in);
    
    // now assign the data matrix
    // double* a = &v[0];
    // maxcol = pos_ind.size();
    maxcol = pos_ptr;
    // cout << "maxcol = " << maxcol << endl;
    Coord_Mat_double Q(read_counter,maxcol,ptr, &qual[0], &read_j[0], &site_i[0], 0); 
    ql_ = Q;
    Coord_Mat_double C(read_counter,maxcol,ptr, &nuc[0], &read_j[0], &site_i[0], 0); 
    seq_ = C;    
    dim_[0]=C.dim(0);
    dim_[1]=C.dim(1);
    nz_ = C.NumNonzeros();
    read_name_ = read_name;
    pos_ = poss;
    ref_=refs;
    alt_=alts;
	return nz_;
}

int HL::read_bam_vcf(const char *bamfile, const char *vcfmarkerfile, const char *hlfile, const char *readfile, const char *sitefile)
// directly write reads' HL from BAM file
// over predefined markers as defined by VCF-formatted file
// VCF file is sorted
// only VCF chr pos ref alt are used
// only expect SNPs
// this is central file for bam_reader
// could be made more efficient by writing out files on the fly
{
	tmpstruct_t tmp;
	tmp.beg = 0; tmp.end = 0x7fffffff;
	tmp.in = samopen(bamfile, "rb", 0);
    tmp.hlfile = fopen( hlfile, "w");
    tmp.readfile = fopen(readfile, "w");
    tmp.sitefile = fopen(sitefile, "w");
    
    
    vector<string> refs;
    vector<string> alts;
    vector<int> poss;
    
    if (tmp.in == 0) {
		fprintf(stderr, "vcf Fail to open BAM file %s\n", bamfile);
		return 1;
	}
    int ref;
    int beg;
    char *ref_allele = new char[1000];
    char *alt_allele = new char[1000];
    bam_index_t *idx;
    bam_plbuf_t *buf;
    idx = bam_index_load(bamfile); // load BAM index
    if (idx == 0) {
        fprintf(stderr, "BAM indexing file is not available.\n");
        return 1;
    }

	char *p = new char [10000]; // Kui Zhang Modified
    
    ifstream myfile;
    myfile.open(vcfmarkerfile);
    string line;
    if (myfile.is_open())
    {
        while ( myfile.good() )
        {
            getline (myfile,line);
            if (line[0] == '#')
                continue;
            if ( line.size()==0 )
                break;
            // cout << "read line: " << line << endl;
            // char *p=new char[line.size()]; // Kui Zhang Modified
            strcpy(p,line.c_str());
            vcf_parse_line(tmp.in->header, p, &ref, &beg, ref_allele, alt_allele); // parse a variant
            tmp.beg = beg-1;
            tmp.end = beg;
            // cout << "ref=" << ref << "; beg=" << tmp.beg << endl;
            if (ref < 0) {
                fprintf(stderr, "Invalid region in file %s\n", vcfmarkerfile);
                return 1;
            }
            
            buf = bam_plbuf_init(pileup_func_nosave, &tmp); // initialize pileup
            bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
            bam_plbuf_push(0, buf); // finalize pileup
        }
        myfile.close();
    }
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);
	
	samclose(tmp.in);
	fclose(tmp.hlfile);
	fclose(tmp.readfile);
	fclose(tmp.sitefile);

    /*
    read_name_ = read_name;
    pos_ = pos_name;
    ref_=refs;
    alt_=alts;
     */
    return 0;
}


int HL::seq_to_s01()
// use ref_ and alt_ to translate seq_ to s01_
{
    // cout << "use ref_ and alt_ to translate seq_ to s01_" <<endl;
    double *s01_val = new double[nz_];
    double *q01_val = new double[nz_];
    int nz = 0;
    int *r = new int[nz_];
    int *c = new int[nz_];
    int col;
    CompRow_Mat_double A=seq_;
    CompRow_Mat_double B=ql_;
    
    for (int i = 0; i < seq_.dim(0) ; i++)
    {
        // char c = ref_[i].c_str()[0];
        for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++) {
            col = A.col_ind(j);
            if ((int)(A.val(j)) == ref_[col].c_str()[0]) {
                s01_val[nz] = 0;
                q01_val[nz] = B.val(j);
            }
            else if ((int)A.val(j) == alt_[col].c_str()[0]) {
                s01_val[nz] = 1;
                q01_val[nz] = B.val(j);
            }
            else {
                // s01_val[nz] = -1;
                // q01_val[nz] = B.val(j);
                continue; // so we ignore non-ref-non-alt bases for now
            }
            r[nz]       = i;
            c[nz]       = col;
            nz ++;
        }
    }
    
    s01_ = Coord_Mat_double(seq_.dim(0),seq_.dim(1),nz,s01_val,r,c);
    q01_ = Coord_Mat_double(seq_.dim(0),seq_.dim(1),nz,q01_val,r,c);
    return nz;
}

int HL::parse_bam(const char *bamfile, const char *markerfile)
// construct HL from BAM file over predefined markers
// bamfile must be sorted and indexed
// markersfile is also sorted
{
	tmpstruct_t tmp;
	tmp.beg = 0; tmp.end = 0x7fffffff;
	tmp.in = samopen(bamfile, "rb", 0);
    
    if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", bamfile);
		return 1;
	}
    int ref;
    bam_index_t *idx;
    bam_plbuf_t *buf;
    idx = bam_index_load(bamfile); // load BAM index
    if (idx == 0) {
        fprintf(stderr, "BAM indexing file is not available.\n");
        return 1;
    }
    
	char *p = new char [10000]; // Kui Zhang Modified
	 
    ifstream myfile; 
    myfile.open(markerfile);
    string line;
    if (myfile.is_open())
    {
        while ( myfile.good() )
        {
            getline (myfile,line);
            if ( line.size()==0 )
                break;
            // cout << "read line: " << line << endl;
            // char *p=new char[line.size()]; // Kui Zhang Modified
            strcpy(p,line.c_str());
            bam_parse_region(tmp.in->header, p, &ref,
                             &tmp.beg, &tmp.end); // parse the region
            if (ref < 0) {
                fprintf(stderr, "Invalid region in file %s\n", markerfile);
                return 1;
            }
            buf = bam_plbuf_init(pileup_func_novcf, &tmp); // initialize pileup
            bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
            bam_plbuf_push(0, buf); // finalize pileup
        }
        myfile.close();
    }
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);
	
	samclose(tmp.in);
    
    // now assign the data matrix
    Coord_Mat_double Q(read_counter,maxcol,ptr, &qual[0], &read_j[0], &site_i[0], 0); 
    ql_ = Q;
    Coord_Mat_double C(read_counter,maxcol,ptr, &qual[0], &read_j[0], &site_i[0], 0); 
    seq_ = C;    
    dim_[0]=C.dim(0);
    dim_[1]=C.dim(1);
    nz_ = C.NumNonzeros();
    read_name_ = read_name;
    pos_ = pos_name;
	return nz_;
}

void HL::write_reads(FILE *out_file)
{
    // write out read index 
    // cout << "write out read index" <<endl;
    fprintf(out_file, "#read index (1-based)\n");
    fprintf(out_file, "#num_of_reads: %d\n", read_name_.size());
    for (int i=0; i < read_name_.size(); i++) {
        fprintf(out_file, "%d\t%s\n", i+1, read_name_[i].c_str());
    }
    fprintf(out_file, "\n");    
}

void HL::write_sites(FILE *out_file, const char *style)
{
    // write out site index
    int bit = 0;
    if (strcmp(strndup(style+strlen(style)-2,2), "01")==0)
        bit = 1;
    
    cout << "write out site index" <<endl;
    fprintf(out_file, "#site index (1-based)\n");
    fprintf(out_file, "#num_of_sites: %d\n", pos_.size());
    // cout << "bit is " << bit << "; style=" << style << endl;
    for (int i=0; i < pos_.size(); i++) {
        // cout << "outputing site " << i << endl;
        if (bit) 
            fprintf(out_file, "%d\t%d\t%s\t%s\n", i+1, pos_[i]+1, ref_[i].c_str(), alt_[i].c_str());            
        else 
            fprintf(out_file, "%d\t%d\n", i+1, pos_[i]+1);
    }
    fprintf(out_file, "\n");
    // cout << "done writeout sizes" << endl;
}

void HL::write_matrix(FILE *out_file, const char *style)
{
    int flag = 0;
    int M = seq_.dim(0);
    int N = seq_.dim(1);
    int rowp1, colp1;
    // write out matrix
    cout << "write out matrix" <<endl;
    fprintf(out_file, "#the matrix (sorted by %s)\n", style);
    fprintf(out_file, "#rid\tsid\tnuc\tqual\n");
    
    if (strcmp(style, "site")==0) {
        CompCol_Mat_double A=seq_;
        CompCol_Mat_double B=ql_;

        //  Loop through columns
        for (int j = 0; j < N ; j++)
            for (int i=A.col_ptr(j);i<A.col_ptr(j+1);i++)
            {
                rowp1 = A.row_ind(i);
                colp1 = j;
                // cout << "i=" << i << "; j=" << j << "; val=" << A.val(i) << endl;
                if ( rowp1 == M-1 && colp1 == N-1 ) flag = 1;
                fprintf(out_file,"%d\t%d\t%c\t%d\n", rowp1+1, colp1, int(A.val(i)), int(B(rowp1,colp1  )));
            }
    }    
    
    if (strcmp(style, "read")==0) {
        CompRow_Mat_double A=seq_;
        CompRow_Mat_double B=ql_;
        
        //  Loop through rows
        for (int i = 0; i < M ; i++)
            for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++)
            {
                colp1 = A.col_ind(j);
                rowp1 = i;
                // cout << "i=" << i << "; j=" << j << "; val=" << A.val(i) << endl;
                if ( rowp1 == M-1 && colp1 == N-1 ) flag = 1;
                fprintf(out_file,"%d\t%d\t%c\t%d\n", rowp1+1, colp1, int(A.val(j)), int(B(rowp1,colp1  )));
            }
    }

    if (strcmp(style, "site01")==0) {
        CompCol_Mat_double A=s01_;
        CompCol_Mat_double B=ql_;
        
        //  Loop through columns
        for (int j = 0; j < N ; j++)
            for (int i=A.col_ptr(j);i<A.col_ptr(j+1);i++)
            {
                rowp1 = A.row_ind(i);
                colp1 = j;
                // cout << "i=" << i << "; j=" << j << "; val=" << A.val(i) << endl;
                if ( rowp1 == M-1 && colp1 == N-1 ) flag = 1;
                fprintf(out_file,"%d\t%d\t%d\t%d\n", rowp1+1, colp1, int(A.val(i)), int(B(rowp1,colp1  )));
            }
    }    
    
    if (strcmp(style, "read01")==0) {
        CompRow_Mat_double A=s01_;
        CompRow_Mat_double B=ql_;
        // CompRow_Mat_double C=seq_;
        
        //  Loop through rows
        for (int i = 0; i < M ; i++)
            for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++)
            {
                colp1 = A.col_ind(j);
                rowp1 = i;
                // cout << "i=" << i << "; j=" << j << "; val=" << A.val(i) << endl;
                if ( rowp1 == M-1 && colp1 == N-1 ) flag = 1;
                fprintf(out_file,"%d\t%d\t%d\t%c\t%d\n", rowp1+1, colp1+1, int(A.val(j)), int(seq_(rowp1,colp1)), int(B(rowp1,colp1  )));
            }
    }

    if (flag == 0) 
        fprintf(out_file, "%d\t%d\t%d\t%d\n", M,N, int(seq_(M-1,N-1)), int(ql_(M-1,N-1)));
}

void HL::write_hl_by_read(const int person_id, const char *outfilename, const int use_quality)
{
    FILE *out_file;
    out_file = fopen( outfilename, "w");
    
    // now write it out
    CompRow_Mat_double A=s01_;
    CompRow_Mat_double B=q01_;
    
    //  Loop through rows
    for (int i = 0; i < A.dim(0) ; i++)
    {
        fprintf(out_file,"%d\t%d\t%d\t%d\t%d", person_id, person_id, 0, 0, 1);
        fprintf(out_file,"\t%d\t%d", i+1, A.row_ptr(i+1)-A.row_ptr(i));
        
        for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++)
            fprintf(out_file,"\t%d", A.col_ind(j)+1);
        for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++)
            fprintf(out_file,"\t%d", int(A.val(j)));
        
        if (use_quality)
            for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++)
                fprintf(out_file,"\t%d", int(B.val(j)));
        
        fprintf(out_file,"\n");
    }
    
    fclose(out_file);
}

void HL::write_counts(const int person_id, const char *outfilename)
{
    FILE *out_file;
    out_file = fopen( outfilename, "w");
    
    // now write it out
    CompCol_Mat_double A=s01_;
    // CompCol_Mat_double B=ql_;
    
    int count01[2];
    //  Loop through columns
    fprintf(out_file,"%d %d %d %d %d", person_id, person_id, 0, 0, 1);
    for (int i = 0; i < A.dim(1) ; i++)
    {
        count01[0]=count01[1]=0;
        for (int j=A.col_ptr(i);j<A.col_ptr(i+1);j++)
            if (int(A.val(j))>=0)
                count01[int(A.val(j))]++;
        fprintf(out_file," %d %d", count01[0], count01[1]);        
    }
    fprintf(out_file,"\n");
    
    fclose(out_file);
}


void HL::write_jumps(const int person_id, const char *of1, const char *of2)
{
    FILE *out_file;
    
    // now write it out
    CompRow_Mat_double A=s01_;
    int b;
    int oddness = 0;
    
    //  simple way of initialize 3D array
    // from http://www.cplusplus.com/forum/articles/7459/
    vector<vector<vector<int> > > count;
    count.resize(2);
    for (int i=0; i<2; i++) {
        count[i].resize(A.dim(1)+1);
    
        for (int j=0; j<=A.dim(1); j++)
            count[i][j].resize(5);
    }
    
    // intialization
    for (int i=0; i<=A.dim(1); i++)
        for (int j=0; j<5; j++)
            for (int k=0; k<2; k++)
            count[k][i][j]=0;
    
    //  Loop through reads (rows)
    for (int i = 0; i < A.dim(0) ; i++)
    {
        int j = A.row_ptr(i);
        if (A.col_ind(j+1)>=A.dim(1)) continue;
        if (A.val(j)<0) continue;
        if (A.val(j+1)<0) continue;
        
        // count num of adjacent pairs
        int pair_stretch_start[100]; // supposedly the max number of pairs possible
        int pair_stretch_length[100]; 
        int pair_stretch_count = -1;
        int in_pair_stretch = 0;
        for (int j=A.row_ptr(i);j<A.row_ptr(i+1)-1;j++) {
            if ((A.col_ind(j+1)-A.col_ind(j)==1) && (!in_pair_stretch) ) {
                pair_stretch_count++;
                pair_stretch_start[pair_stretch_count] = j;
                pair_stretch_length[pair_stretch_count] = 1;
                in_pair_stretch = 1;
                continue;
            }
            if ((A.col_ind(j+1)-A.col_ind(j)==1) && (in_pair_stretch) ) {
                pair_stretch_length[pair_stretch_count] ++;
                // cout << j << " " << "in pair stretch " << pair_stretch_count << " len: " << pair_stretch_length[pair_stretch_count] << endl;
		continue;
            }
            else
            {
                in_pair_stretch = 0;
            }
        }
        pair_stretch_count++;
        
        for (int ps=0; ps<pair_stretch_count; ps++)
        {
            if (pair_stretch_length[ps] == 1) {
                j = pair_stretch_start[ps];
                b = int(A.val(j)*2+A.val(j+1));
                count[0][A.col_ind(j)][b] ++;
                count[0][A.col_ind(j)][4] ++; // contains the total base count
                count[1][A.col_ind(j)][b] ++;
                count[1][A.col_ind(j)][4] ++; // contains the total base count
                continue;
            } else {
                // pairs_count >=2
                oddness = 0;
                j = pair_stretch_start[ps];
                // cout << j << " length: " << pair_stretch_length[ps] << endl;
		for (int pi=0; pi < pair_stretch_length[ps]; pi++) {
                    b = int(A.val(j)*2+A.val(j+1));
                    count[oddness][A.col_ind(j)][b] ++;
                    count[oddness][A.col_ind(j)][4] ++; // contains the total base count
                    oddness = 1-oddness; // flip oddness
                    j ++;
                }
            }
        }
    }
    for (oddness=0; oddness < 2; oddness ++) {
        if (oddness == 0) out_file = fopen( of1, "w");
            else out_file = fopen( of2, "w");
        
        for (int i = 0; i < A.dim(1)-1 ; i++)
        {
            if (count[oddness][i][4]==0) continue;
            fprintf(out_file,"%d\t%d\t%d\t%d\t%d", person_id, person_id, 0, 0, 1);
            fprintf(out_file,"\t%d\t%d", i+1, i+2); // 1-based
            for (b=0; b<4; b++) {
                fprintf(out_file,"\t%d", count[oddness][i][b]);
            }
            fprintf(out_file,"\n");
        }
        fclose(out_file);
    }
    
}



void HL::write(const char *filename, const char *style)
{
    char my_style[100];
    strcpy(my_style, style);
    if (s01_.NumNonzeros()) strcat(my_style, "01");
    FILE *out_file;
    out_file = fopen( filename, "w");
    
    write_reads(out_file);
    write_sites(out_file, my_style);
    // write_matrix(out_file, my_style);
    fclose(out_file);
    
    return;
}

int		HL::get_single_reads(vector<int>& rs) 
// return the number of single-site reads
{
    // int * rs = new int[1000]; // need allocate memory
    int subset_ind = 0;
    rs.clear();
    rs.resize(1000);
    CompRow_Mat_double A=seq_;
    
    for (int i = 0; i < seq_.dim(0) ; i++)
        if (A.row_ptr(i+1)-A.row_ptr(i) == 1)
            rs[subset_ind++] = i;
    rs.resize(subset_ind);
    return subset_ind;
}

int		HL::get_long_reads(vector<int>& rs) 
// return the number of single-site reads
{
    // int * rs = new int[1000]; // need allocate memory
    // int * rs = new int[1000]; // need allocate memory
    int subset_ind = 0;
    rs.resize(1000);
    CompRow_Mat_double A=seq_;
    
    for (int i = 0; i < seq_.dim(0) ; i++)
        if (A.row_ptr(i+1)-A.row_ptr(i) > 2)
            rs[subset_ind++] = i;
    rs.resize(subset_ind);
    return subset_ind;
}

int		HL::get_double_reads(vector<int>&  rs) 
// return the number of double-site reads
{
    // read_subset = new int[1000]; // need allocate memory
    // int * rs = new int[1000]; // need allocate memory
    int subset_ind = 0;
    rs.resize(1000);
    CompRow_Mat_double A=seq_;
    
    for (int i = 0; i < seq_.dim(0) ; i++)
        if (A.row_ptr(i+1)-A.row_ptr(i) == 2)
            rs[subset_ind++] = i;
    rs.resize(subset_ind);
    return subset_ind;
}

HL HL::subset_by_sites(vector<int> pos_ind)
{
    HL h;
    return h;
}

HL HL::subset_by_reads(int *read_ind, int M)
{
    HL h;
    CompRow_Mat_double A=seq_;
    CompRow_Mat_double B=ql_;
    CompRow_Mat_double C=s01_;
    CompRow_Mat_double D=q01_;
    
    // prepare seq_val_, ql_val_, nz_, M, N, r, c, base=0
    int N=A.dim(1);
    
    int read_ind_ptr = 0;
    int nz = 0;
    double *seq_val = new double[nz_];
    double *ql_val = new double[nz_];
    double *s01_val = new double[nz_];
    double *q01_val = new double[nz_];
    int *r = new int[nz_];
    int *c = new int[nz_];
    vector<string> read_name;
    read_name.reserve(10000);
    
    cout << "in int*" << endl;
    for (int i = 0; i < seq_.dim(0) ; i++)
    {
        if (i==read_ind[read_ind_ptr]) {
            read_name[read_ind_ptr] = read_name_[i];
            read_ind_ptr++;
            if (read_ind_ptr % 10000 ==0)
                read_name.reserve(read_ind_ptr+10000);
        }
        else
            continue;
        
        for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++) {
            seq_val[nz] = A.val(j);
            ql_val[nz]  = B.val(j);
            s01_val[nz] = C.val(j);
            q01_val[nz] = D.val(j);
            r[nz]       = read_ind_ptr-1;
            c[nz]       = A.col_ind(j);
            nz ++;
        }
    }
    
    h.seq_ = Coord_Mat_double(M,N,nz,seq_val,r,c);
    h.ql_  = Coord_Mat_double(M,N,nz,ql_val ,r,c);
    h.s01_ = Coord_Mat_double(M,N,nz,s01_val,r,c);
    h.q01_ = Coord_Mat_double(M,N,nz,q01_val,r,c);
    h.nz_ = nz;

    h.dim_[0]=h.ql_.dim(0);
    h.dim_[1]=h.ql_.dim(1);
//    std::copy(read_name, read_name+read_ind_ptr, std::back_inserter(h.read_name_));
    h.read_name_ = read_name;
    
    h.pos_ = pos_;
    h.ref_ = ref_;
    h.alt_ = alt_;
    /*int *pos = new int[N];
    string *ref = new string[N];
    string *alt = new string[N];
    for (int j=0; j< N; j++)
    {
        pos[j] = pos_[j];
        ref[j] = ref_[j];
        alt[j] = alt_[j];
    }
    
    h.pos_ = pos;
    h.ref_ = ref;
    h.alt_ = alt;
     */
    return h;
}

HL HL::subset_by_reads(vector<int>& read_ind, int M)
{
    HL h;
    CompRow_Mat_double A=seq_;
    CompRow_Mat_double B=ql_;
    CompRow_Mat_double C=s01_;
    CompRow_Mat_double D=q01_;
    
    // prepare seq_val_, ql_val_, nz_, M, N, r, c, base=0
    int N=A.dim(1);
    
    int read_ind_ptr = 0;
    int nz = 0;
    double *seq_val = new double[nz_];
    double *ql_val = new double[nz_];
    double *s01_val = new double[nz_];
    double *q01_val = new double[nz_];
    int *r = new int[nz_];
    int *c = new int[nz_];
    vector<string> read_name;
    read_name.reserve(10000);
    
    cout << "in vector. size of read_ind" << read_ind.size() << endl;
    if (read_ind.size()) {
    for (int i = 0; i < seq_.dim(0) ; i++)
    {
        if (i==read_ind[read_ind_ptr]) {
            read_name[read_ind_ptr] = read_name_[i];
            read_ind_ptr++;
            if (read_ind_ptr % 10000 ==0)
                read_name.reserve(read_ind_ptr+10000);
        }
        else
            continue;
        
        for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++) {
            seq_val[nz] = A.val(j);
            ql_val[nz]  = B.val(j);
            s01_val[nz] = C.val(j);
            q01_val[nz] = D.val(j);
            r[nz]       = read_ind_ptr-1;
            c[nz]       = A.col_ind(j);
            nz ++;
        }
    }
    }
    cout << "nz=" << nz << endl;
    
    h.seq_ = Coord_Mat_double(M,N,nz,seq_val,r,c);
    h.ql_  = Coord_Mat_double(M,N,nz,ql_val ,r,c);
    h.s01_ = Coord_Mat_double(M,N,nz,s01_val,r,c);
    h.q01_ = Coord_Mat_double(M,N,nz,q01_val,r,c);
    h.nz_ = nz;
    
    h.dim_[0]=h.ql_.dim(0);
    h.dim_[1]=h.ql_.dim(1);
    h.read_name_ = read_name;
    
    h.pos_ = pos_;
    h.ref_ = ref_;
    h.alt_ = alt_;
    
/*    int *pos = new int[N];
    string *ref = new string[N];
    string *alt = new string[N];
    for (int j=0; j< N; j++)
    {
        pos[j] = pos_[j];
        ref[j] = ref_[j];
        alt[j] = alt_[j];
    }
    
    h.pos_ = pos;
    h.ref_ = ref;
    h.alt_ = alt;
  */  
    return h;
}

HL HL::subset(vector<int> pos_ind, vector<int> reads_ind)
{
    HL h;
    return h;
}


double** HL::calculate_GL()
{
    double **GL = new double *[3];
    
    for (int i=0; i<3; i++) {
        GL[i] = new double[dim(1)];
        for (int j=0; j<dim(1); j++)
            GL[i][j]=0;
    }
    
    CompCol_Mat_double A=seq_;
    CompCol_Mat_double B=ql_;

    //  Loop through columns
    for (int j = 0; j < dim(1) ; j++) {
        // first get the major and minor alleles
        for (int i=A.col_ptr(j);i<A.col_ptr(j+1);i++)
        {
            int rowp1 = A.row_ind(i);
            int colp1 = j;
            // cout << "i=" << i << "; j=" << j << "; val=" << A.val(i) << endl;
        }
        // then calculate GLs
        for (int i=A.col_ptr(j);i<A.col_ptr(j+1);i++)
        {
            int rowp1 = A.row_ind(i);
            int colp1 = j;
            // cout << "i=" << i << "; j=" << j << "; val=" << A.val(i) << endl;
        }
    }
    return GL;
}

vector<int> HL::merge_read_sets(const vector<int>& read_subset1, const vector<int>& read_subset2)
// merge two sorted indices
{
    vector<int> rs(read_subset1.size()+read_subset2.size());
       
    merge(read_subset1.begin(), read_subset1.end(), read_subset2.begin(), read_subset2.end(), rs.begin());
    return rs;
}

int HL::split(HL& h, int oddness)	
{
    CompRow_Mat_double A=seq_;
    CompRow_Mat_double B=ql_;
    CompRow_Mat_double C=s01_;
    CompRow_Mat_double D=q01_;
    
    // prepare seq_val_, ql_val_, nz_, M, N, r, c, base=0
    int N=A.dim(1);
    
    int read_ind_ptr = 0;
    int nz = 0;
    int j, k;
    char k_string[100];
    double *seq_val = new double[nz_];
    double *ql_val = new double[nz_];
    double *s01_val = new double[nz_];
    double *q01_val = new double[nz_];
    int *r = new int[nz_];
    int *c = new int[nz_];
    vector<string> read_name;
    read_name.reserve(10000);
    
    cout << "in split. oddness=" << oddness << endl;
    if (read_ind.size()) { //// problem read_ind undefined!
        for (int i = 0; i < A.dim(0) ; i++)
        {
            if (A.row_ptr(i+1)-A.row_ptr(i)<=2)
            {
                for (j=A.row_ptr(i);j<A.row_ptr(i+1);j++) {
                    read_name[read_ind_ptr] = read_name_[i];
                    seq_val[nz] = A.val(j);
                    ql_val[nz]  = B.val(j);
                    s01_val[nz] = C.val(j);
                    q01_val[nz] = D.val(j);
                    r[nz]       = read_ind_ptr;
                    c[nz]       = A.col_ind(j);
                    nz ++;
                    read_ind_ptr ++;
                    if (read_ind_ptr % 10000 ==0)
                        read_name.reserve(read_ind_ptr+10000);
                }
                continue;
            }
            
            // now split is needed
            j=A.row_ptr(i);
            k = 0; // the index of split segment
            if (oddness) {
                sprintf(k_string, "%d", k);
                read_name[read_ind_ptr] = read_name_[i]+string(k_string);
                seq_val[nz] = A.val(j);
                ql_val[nz]  = B.val(j);
                s01_val[nz] = C.val(j);
                q01_val[nz] = D.val(j);
                r[nz]       = read_ind_ptr;
                c[nz]       = A.col_ind(j);
                nz ++;
                read_ind_ptr ++;                
            }
            while (j+2<A.row_ptr(i+1))
            {
                for (int m=0;m<2;m++) {
                    k++;
                    j++;
                    sprintf(k_string, "%d", k);
                    read_name[read_ind_ptr] = read_name_[i]+string(k_string);
                    seq_val[nz] = A.val(j);
                    ql_val[nz]  = B.val(j);
                    s01_val[nz] = C.val(j);
                    q01_val[nz] = D.val(j);
                    r[nz]       = read_ind_ptr;
                    c[nz]       = A.col_ind(j);
                    nz ++;
                }
                read_ind_ptr ++;                                
            }
            
            if (j+2==A.row_ptr(i+1))
            {
                k ++;
                j ++;
                sprintf(k_string, "%d", k);
                read_name[read_ind_ptr] = read_name_[i]+string(k_string);
                seq_val[nz] = A.val(j);
                ql_val[nz]  = B.val(j);
                s01_val[nz] = C.val(j);
                q01_val[nz] = D.val(j);
                r[nz]       = read_ind_ptr;
                c[nz]       = A.col_ind(j);
                nz ++;
                read_ind_ptr ++;                                
            }
        }
    }
    cout << "nz=" << nz << "nz_=" << nz_ << endl;
    
    h.seq_ = Coord_Mat_double(read_ind_ptr,N,nz,seq_val,r,c);
    h.ql_  = Coord_Mat_double(read_ind_ptr,N,nz,ql_val ,r,c);
    h.s01_ = Coord_Mat_double(read_ind_ptr,N,nz,s01_val,r,c);
    h.q01_ = Coord_Mat_double(read_ind_ptr,N,nz,q01_val,r,c);
    h.nz_ = nz;
    
    h.dim_[0]=h.ql_.dim(0);
    h.dim_[1]=h.ql_.dim(1);
    h.read_name_ = read_name;
    
    // simply copy sites    
    h.pos_ = pos_;
    h.ref_ = ref_;
    h.alt_ = alt_;
    
    return read_ind_ptr;
    
    // merge back single, double, and split long reads in to even an odd valition,
    
}

