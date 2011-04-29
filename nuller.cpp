#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <iostream>
#include "convert.h"
#include "api/BamMultiReader.h"
#include "fastahack/Fasta.h"
#include "vcflib/Variant.h"

using namespace BamTools;
using namespace std;

void printUsage(int argc, char** argv) {

    cerr << "usage: " << argv[0] << " [options]" << endl
         << endl
         << "options:" << endl
         << "    -h, --help         this dialog" << endl
         << "    -b, --bam FILE     use this BAM as input (any number, or none for stdin)" << endl
         << "    -q, --min-bq QUAL  threshold PHRED scaled base quality" << endl
         << "    -m, --min-mq QUAL  threshold PHRED scaled mapping quality" << endl
         << endl
         << "Outputs a VCF where samples are called as reference where they have data." << endl
         << "The output is intended for use in data availability analysis." << endl
         << endl
         << "author: Erik Garrison <erik.garrison@bc.edu>" << endl;

}

void writeNullVCFHeader(ostream& out, vector<string>& samples, string& referenceName, string& commandline) {

    out << "##fileformat=VCFv4.0" << endl
        << "##source=nuller" << endl
        << "##reference=" << referenceName << endl
        << "##phasing=none" << endl
        << "##commandline=\"" << commandline << "\"" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl
        << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
            out << "\t" << *s;
        }

    out << endl;

}

void writeNullVCFEntry(ostream& out,
        string& sequence,
        int position,
        FastaReference& reference,
        vector<string>& sampleList,
        map<string, int>& coverage) {

    string refbase = reference.getSubSequence(sequence, position, 1);

    int depth = 0;
    for (map<string, int>::iterator c = coverage.begin(); c != coverage.end(); ++c) {
        depth += c->second;
    }

    out << sequence << "\t"
        << position << "\t"
        << "." << "\t"
        << refbase << "\t"
        << refbase << "\t"
        << "0" << "\t"
        << "." << "\t"
        << "NS=" << coverage.size() << ";DP=" << depth << "\t"
        << "GT:DP";

    for (vector<string>::iterator s = sampleList.begin(); s != sampleList.end(); ++s) {
        string& name = *s;
        map<string, int>::iterator f = coverage.find(name);
        if (f != coverage.end()) {
            out << "\t" << "./.:" << f->second;
        } else {
            out << "\t" << ".";
        }
    }

    out << endl;

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;
    int minBaseQuality = 0;
    int minMappingQuality = 0;
    string fastaReferenceFilename;
    FastaReference reference;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"min-bq", required_argument, 0, 'q'},
            {"min-mq", required_argument, 0, 'm'},
            {"fasta", required_argument, 0, 'f'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:q:m:f:",
                         long_options, &option_index);

        if (c == -1)
            break;
 
        switch (c) {

            case '?':
                printUsage(argc, argv);
                exit(1);
                break;

            case 'h':
                printUsage(argc, argv);
                exit(1);
                break;

            case 'b':
                inputFilenames.push_back(optarg);
                break;

            case 'q':
                if (!convert(optarg, minBaseQuality)) {
                    cerr << "could not convert --min-mapq" << endl;
                    exit(1);
                }
                break;

            case 'm':
                if (!convert(optarg, minMappingQuality)) {
                    cerr << "could not convert --min-baseq" << endl;
                    exit(1);
                }
                break;

            case 'f':
                fastaReferenceFilename = optarg;
                break;

            default:
                exit(1);
                break;
        }
    }

    if (inputFilenames.empty()) {
        inputFilenames.push_back("stdin");
    }

    if (fastaReferenceFilename.empty()) {
        cerr << "no fasta reference supplied.  but do we really need it?" << endl;
        exit(1);
    } else {
        reference.open(fastaReferenceFilename);
    }

    BamMultiReader reader;
    if (!reader.Open(inputFilenames)) {
        cerr << "could not open input BAM files" << endl;
        exit(1);
    }

    RefVector references = reader.GetReferenceData();
    SamHeader header = reader.GetHeader();

    map<int, string> refIDtoSequenceName;
    int i = 0;
    for (RefVector::const_iterator r = references.begin(); r != references.end(); ++r) {
        refIDtoSequenceName[i++] = r->RefName;
    }

    vector<string> sampleList;
    map<string, string> samplesToReadGroupIDs;
    for (SamReadGroupConstIterator g = header.ReadGroups.Begin(); g != header.ReadGroups.End(); ++g) {
        samplesToReadGroupIDs[g->Sample] = g->ID;
    }
    
    for (map<string, string>::iterator s = samplesToReadGroupIDs.begin(); s != samplesToReadGroupIDs.end(); ++s) {
        sampleList.push_back(s->first);
    }

    string commandline;
    for (int i = 0; i < argc; ++i) {
        commandline += " " + string(argv[i]);
    }

    // write the VCF header
    writeNullVCFHeader(cout, sampleList, fastaReferenceFilename, commandline);

    map<int, map<string, int> > sampleCoverageMap;
    int currRefID = 0;
    string currSequenceName = refIDtoSequenceName[currRefID];
    int currPosition = 0;
    int lastErasedPosition = 0; // for cleaning up the map

    BamAlignment alignment;
    string readGroup;
    while (reader.GetNextAlignment(alignment)) {

        if (!alignment.IsMapped() || alignment.IsDuplicate() || !alignment.IsPrimaryAlignment())
            continue;

        if (alignment.MapQuality < minMappingQuality)
            continue;

        if (alignment.RefID != currRefID) {
            // we've switched references, so dump the last reference
            for (map<int, map<string, int> >::iterator s = sampleCoverageMap.begin();
                    s != sampleCoverageMap.end(); ++s) {
                int position = s->first;
                map<string, int>& coverageMap = s->second;
                writeNullVCFEntry(cout,
                        currSequenceName,
                        position,
                        reference,
                        sampleList,
                        coverageMap);
            }
            sampleCoverageMap.clear();
            currRefID = alignment.RefID;
            currSequenceName = refIDtoSequenceName[currRefID];
        }

        if (alignment.Position > currPosition) {
            // dump everything up to the current position, as we can't get
            // more coverage information there
            currPosition = alignment.Position;
            for (int i = lastErasedPosition; i < currPosition; ++i) {
                map<int, map<string, int> >::iterator f = sampleCoverageMap.find(i);
                if (f != sampleCoverageMap.end()) {
                    int position = f->first;
                    map<string, int>& coverageMap = f->second;
                    writeNullVCFEntry(cout,
                            currSequenceName,
                            position,
                            reference,
                            sampleList,
                            coverageMap);
                    sampleCoverageMap.erase(i);
                }
            }
            lastErasedPosition = currPosition;
        }

        if (alignment.GetTag("RG", readGroup)) {
            string sampleName = header.ReadGroups[readGroup].Sample;
            // parse the cigar, quality string

            alignment.QueryBases;

            int rp = 0;  // read position, 0-based relative to read
            int sp = alignment.Position;  // sequence position

            for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin();
                c != alignment.CigarData.end(); ++c) {
                unsigned int l = c->Length;
                char t = c->Type;
                if (t == 'M') { // match or mismatch
                    for (int i = 0; i < l; ++i, ++sp) {
                        ++sampleCoverageMap[sp][sampleName];
                    }
                    rp += l;
                } else if (t == 'D') { // deletion
                    sp += l;  // update reference sequence position
                } else if (t == 'I') { // insertion
                    rp += l;
                } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
                    rp += l;
                } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
                } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
                    sp += l;
                }
            }

        }

    }

    reader.Close();

    return 0;

}
