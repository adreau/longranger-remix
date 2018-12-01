
#####################################################################################
# version 3, no modeling of contamination, DO sanity check of the neighbor region ###
#                                                                                 ###
#####################################################################################

import longranger.cnv.rpm_utils as ru
from  shutil import copyfileobj
import os

__MRO__ = """
stage CALL_CNV(
    in  string[]  rPMFiles,
    in  int       maxUnknownGap,
    in  int       maxNormalGap,
    in  int       maxHomoGap,
    in  bool      wgsmode,
    out bed       het_del,
    src py       "stages/cnv/call_cnv",
)
"""

LENGTH_CAP=2000000

def split(args):
    if args.rPMFiles:
        files = [{"rPM":rpm, "__mem_gb": 12.0, '__threads':2} for rpm in args.rPMFiles]
        return {'chunks': files}
    else:
        return {'chunks':[{"rPM":None}]}


def main(args, outs):

    if not args.rPMFiles:
        outs.het_del = None
        return

    rpm_by_chrom = ru.ReadRPM(args.rPM)

    het_del=open(outs.het_del,'w')
    print outs.het_del

    sorted_chrom = sorted(rpm_by_chrom.keys())
    for chrom in sorted_chrom:
        bins = rpm_by_chrom[chrom]
        IDs = sorted(bins.keys())
        ln = len(IDs)
        binIDMin = IDs[0]
        binIDMax = IDs[ln-1]

        idx = -1
        while idx < ln-1:

            idx += 1
            ID = IDs[idx]
            r = bins[ID]
            if r.exonID %10 ==1: continue # off-target region
            hasMissing, whichMissing = r.testMissingHP()
            #print r.chrom, r.start, r.end, r.exonID, ID

            if hasMissing == 1:
                # initialize a het_del test
                id_list = [ID]
                hp_list = [whichMissing]
                lastHet = r # the end of the last bin showing sign of het del
                lastMissing = whichMissing

                rLoc = r
                while idx < ln-1:
                    idx += 1
                    IDCur = IDs[idx]
                    rPre = rLoc
                    rLoc = bins[IDCur]
                    hasMissing, whichMissing = rLoc.testMissingHP()
                    #print hasMissing, whichMissing, lastMissing, r[IDs[idx]][1], r[IDs[idx]][0] - r[IDs[idx-1]][1]

                    # when the region is extended:
                    #    1. next bin is the same type and the gap is smaller than 200 bp
                    #    2. reads on both haps for the next bin, next bin is adjacenct to the previous one, the strech is no larger than 200 bp
                    #    3. "homo" type bin, gap is no larger than 200 bp, the strech can not be largert than 1kb??
                    #    4. NO extension when the missing hap got swtiched, say from hap1 to hap2.
                    #
                    #    a strech starts from the end pos of the last het  lastHet
                    #
                    #    if normal and homo bins alternate, whichever breaks the condition break the region continuation


                    #print rLoc.start, rLoc.end, hasMissing, whichMissing, rLoc.exonID, IDCur

                    gapToPre =  rLoc.start - rPre.end
                    if gapToPre > args.maxUnknownGap:
                        print "unknown gap"
                        #print "breaking due to gapToPre", IDs[idx]
                        break

                    if args.wgsmode:
                        if hasMissing == 1: # het del
                            if whichMissing != lastMissing:
                                break
                        elif hasMissing == 0: # homo
                            if rLoc.end - lastHet.end > args.maxHomoGap:
                                break
                        elif hasMissing == -1: # none
                            if rLoc.end - lastHet.end > args.maxNormalGap:
                                break
                    else:
                        if rLoc.exonID % 10 ==1:  # off-target region
                            if hasMissing ==-1:
                                break
                        if hasMissing == -1 or (hasMissing==1 and whichMissing != lastMissing) or abs(rLoc.exonID- lastHet.exonID) > 10:
                            break

                    if hasMissing ==1 and whichMissing == lastMissing:
                        id_list.append(IDCur)
                        hp_list.append(whichMissing)
                        lastHet = rLoc

                print
                print id_list
                for id in id_list:
                    r_tmp = bins[id]
                    print r_tmp.chrom, r_tmp.start, r_tmp.end, r_tmp.exonID, id

                # trim off if the last several bins with reads on both
                #backward = len(hp_list) -1
                #while hp_list[backward]==-1:
                #    backward -= 1
                #id_list=id_list[:backward+1]
                #hp_list=hp_list[:backward+1]

                mergedBin = ru.MergedBin([bins[id] for id in id_list])
                print "merged"
                print mergedBin.chrom, mergedBin.start, mergedBin.end, mergedBin.WildMol, mergedBin.WildMolMs, mergedBin.WildRead, mergedBin.DelMol

                ###############################
                # check for background sanity #
                ###############################

                ### find the range of bin id where id % 5 ==0
                ### find the bins with id % 10 ==0 and size > 200 bp
                ### merge the info and calculat some statistics

                MinBinID, MaxBinID = id_list[0], id_list[len(id_list)-1] # min and max of the binIDs
                MinBinExonID, MaxBinExonID = bins[MinBinID].exonID, bins[MaxBinID].exonID # min and max of the exon IDs
                MinFlankingID, MaxFlankingID = ru.findFlankingIDRange(MinBinExonID, MaxBinExonID)

                ##### find the list of flanking bins
                flanking_left_id = None
                flanking_right_id = None
                ### find the left most flanking bins
                i = MinBinID -1

                isMinFlankingFound = False
                isMaxFlankingFound = False
                while i >= binIDMin and (not isMinFlankingFound):
                    if i in bins:
                        if bins[i].exonID%10==0 and bins[i].end - bins[i].start>=200:
                            if bins[i].exonID <= MinFlankingID:
                                isMinFlankingFound = True
                            flanking_left_id = i
                    i-=1

                i = MaxBinID
                while i <= binIDMax and (not isMaxFlankingFound):
                    if i in bins:
                        if bins[i].exonID%10==0 and bins[i].end - bins[i].start >=200:
                            if bins[i].exonID >= MaxFlankingID:
                                isMaxFlankingFound = True
                            flanking_right_id = i
                    i+=1

                if flanking_left_id and flanking_right_id:
                    flanking_ids = [flanking_left_id, flanking_right_id]
                elif flanking_left_id and (not flanking_right_id):
                    flanking_ids = [flanking_left_id]
                elif (not flanking_left_id) and flanking_right_id:
                    flanking_ids = [flanking_right_id]
                else:
                    flanking_ids = []

                info = ru.calBackgroundStats(bins, flanking_ids)

                deleted_haplotype = hp_list[0] - 1
                if mergedBin.ttlMol ==0: continue
                if mergedBin.hp1Mol != mergedBin.hp1MolMs and mergedBin.hp2Mol != mergedBin.hp2MolMs: continue
                if mergedBin.hp1Mol==0 and mergedBin.hp2Mol==0: continue
                if mergedBin.HetDelPVal <= 1e-4 and mergedBin.WildRead >= 10 and mergedBin.end - mergedBin.start <= LENGTH_CAP:
                    het_del.write("%s\t%d\t%d\t1.0\t%.3e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" %
                            (mergedBin.chrom, mergedBin.start, mergedBin.end, mergedBin.HetDelPVal, \
                             deleted_haplotype, mergedBin.WildRead,mergedBin.size,\
                             mergedBin.ttlMol, mergedBin.ttlMolMs, mergedBin.WildMol, mergedBin.WildMolMs,\
                                 mergedBin.DelMol, mergedBin.DelMolMs, info))

                ### reset idx
                lastID = id_list[len(id_list)-1]
                i = idx;
                while i>=0:
                    if IDs[i]==lastID: break
                    i-=1
                idx = i

    het_del.close()


def join(args, outs, chunk_defs, chunk_outs):
    hetDelFiles = [out.het_del for out in chunk_outs if out.het_del and os.path.exists(out.het_del) ]
    with open(outs.het_del,"w") as fout:
        for infile in hetDelFiles:
            with open(infile, "rb") as fin:
                copyfileobj(fin, fout)
