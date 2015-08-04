//
//  bam_reader.cpp
//  HL
//
//  Created by Degui Zhi on 10/30/12.
//
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
        cerr << "Usage: bam_reader bamfile VCF_markerfile outfile-prefix [person_id]" << endl;
        exit(-1);
    }
    char outprefix[1000], logfile[1000], out_count[1000], out_j1[1000], out_j2[1000], out_read[1000], readlogfile[1000], sitelogfile[1000];
    
    // cout << argv[3] << endl;
    if (argc >= 3)
        strcpy(outprefix, argv[3]);
    else
        strcpy(outprefix, "out");
    sprintf(logfile, "%s.log", outprefix);
    sprintf(readlogfile, "%s.read.log", outprefix);
    sprintf(sitelogfile, "%s.site.log", outprefix);
    sprintf(out_count, "%s.count.txt", outprefix);
    sprintf(out_j1, "%s.jump1.txt", outprefix);
    sprintf(out_j2, "%s.jump2.txt", outprefix);
    sprintf(out_read, "%s.read.txt", outprefix);
    
    
    int person_id = 1;
    if (argc >= 4) person_id = atoi(argv[4]);
    
    HL h;
    int a = h.read_bam_vcf(argv[1], argv[2], out_read, readlogfile, sitelogfile); //
    // cout << "loaded " << a << " entries\n" <<endl;
    h.write(logfile, "site");
   
	return 0;
}
