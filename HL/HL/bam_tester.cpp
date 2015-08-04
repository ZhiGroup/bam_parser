//
//  bam_parser.cpp
//  HL
//
//  Created by Degui Zhi on 8/2/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream.h>
#include <fstream.h>
#include <string.h>
#include "HL.h"
#include "sam.h"
#include "iotext_double.h"

using namespace std;


int main(int argc, char *argv[])
{
    if (argc < 2) {
        cerr << "Usage: bam_tester bamfile VCF_markerfile" << endl;
        exit(-1);
    }
    
    HL h;
    int a = h.parse_bam_vcf(argv[1],argv[2]);
    h.seq_to_s01();
    // int a = h.parse_bam(argv[1],argv[2]);
    cout << "loaded " << a << " entries\n" <<endl;
    // writetxtfile_mat("seq.r.mat.txt",h.seq_Row_);
    // writetxtfile_mat("seq.c.mat.txt",h.seq_Col_);
    //writetxtfile_mat("qual.mat.txt",h.seq_Row_);
    h.write("seq.txt", "read");
    
    // int *double_reads;
    vector<int> single_reads, double_reads;
    // double_reads = new int[1000];
    cout << "number of reads: " << h.dim(0) << endl;
    int k = h.get_single_reads(single_reads);
    cout << "number of single reads: " << k  << endl;
        
    HL h_sub = h.subset_by_reads(single_reads, k);
    cout << "number of single reads in h_s: " << h_sub.get_single_reads(single_reads) << endl;
    cout << "number of double reads in h_s: " << h_sub.get_double_reads(double_reads) << endl;
    h_sub.write("seq_single.txt", "read");

    k = h.get_double_reads(double_reads);
    cout << "number of double reads: " << k << endl;
    HL h_dub = h.subset_by_reads(double_reads, k);
    cout << "number of single reads in h_d: " << h_dub.get_single_reads(single_reads) << endl;
    cout << "number of double reads in h_d: " << h_dub.get_double_reads(double_reads) << endl;
    h_dub.write("seq_double.txt", "read");

	return 0;
}
