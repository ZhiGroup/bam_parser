//
//  bam_parser.cpp
//  HL
//
//  Created by Degui Zhi on 8/2/12.
//  Copyright (c) 2012-2013 __DeguiZhi__. All rights reserved.
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
        cerr << "Usage: bam_parser bamfile VCF_markerfile outfile-prefix [person_id] [use_quality] [use_log]" << endl;
        exit(-1);
    }
    char outprefix[1000], logfile[1000], out_count[1000], out_j1[1000], out_j2[1000], out_read[1000];

    // cout << argc << endl;
    if (argc > 3)
        strcpy(outprefix, argv[3]);
    else
        strcpy(outprefix, "out");
    sprintf(logfile, "%s.log", outprefix);
    sprintf(out_count, "%s.count.txt", outprefix);
    sprintf(out_j1, "%s.jump1.txt", outprefix);
    sprintf(out_j2, "%s.jump2.txt", outprefix);
    sprintf(out_read, "%s.read.txt", outprefix);
    
    
	cout << "bam_parser started...\n";
    int person_id = 1;
    if (argc > 4) person_id = atoi(argv[4]);

    int use_quality = 0;
    if (argc > 5) use_quality = atoi(argv[5]);

    int use_log = 0;
    if (argc > 6) use_log = atoi(argv[6]);
    
    HL h;
    int a = h.parse_bam_vcf(argv[1],argv[2]);
    int b = h.seq_to_s01();
    cout << "loaded " << a << " entries; " << a-b <<  " of them are non-ref-non-alt\n" <<endl;
    
    h.write_hl_by_read(person_id, out_read, use_quality);
    if (use_quality) {
        cout << "done writing read files with quality values" <<endl;
    } else {
        h.write_counts(person_id, out_count);
        h.write_jumps(person_id, out_j1, out_j2);
        cout << "done writing count, jump, and read files" <<endl;
    }
    
    if (use_log) h.write(logfile, "read");

	return 0;
}
