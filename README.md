# Long Ranger (modified for ReMIX)

[Long Ranger](https://github.com/10XGenomics/longranger) pipeline modified for the requirements of [ReMIX](https://github.com/adreau/ReMIX).
A new sub-pipeline called Long Ranger ReportMolecules was added.
This sub-pipeline is based on two parts of the Long Ranger Whole Genome Phasing and SV Calling pipeline (Long Ranger wgs): the computational reconstruction of the molecules and the report of the molecule information in the info field of the vcf file.
This sub-pipeline requires in the input:
- the bam file sorted by genomic positions
- the bam file sorted by barcodes
- the vcf file
- the reference sequence


## License

The code is licensed under the [GNU Affero General Public License version 3](https://support.10xgenomics.com/docs/license).
