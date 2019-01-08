# Long Ranger (modified for ReMIX)

[Long Ranger](https://github.com/10XGenomics/longranger) pipeline modified for the requirements of [ReMIX](https://github.com/adreau/ReMIX).
A new sub-pipeline called Long Ranger ReportMolecules was added.
This sub-pipeline is based on two parts of the Long Ranger Whole Genome Phasing and SV Calling pipeline (Long Ranger wgs): the computational reconstruction of the molecules and the report of the molecule information in the info field of the vcf file.
This sub-pipeline requires in the input:
- the bam file sorted by genomic positions
- the bam file sorted by barcodes
- the vcf file
- the reference sequence

# Usage

- Download and unpack the [Long Ranger](https://support.10xgenomics.com/genome-exome/software/downloads/latest) archive.

- Clone or download the *longranger-remix* repository into *longranger-[version]/longranger-cs* folder
```sh
cd longranger-2.2.2/longranger-cs
git clone https://github.com/adreau/longranger-remix.git
```

- Create a symbolic link *longranger-[version]/longranger-remix* to *longranger-[version]/longranger-cs/longranger-remix/bin/longranger*
```sh
cd ..
ln -s longranger-cs/longranger-remix/bin/longranger longranger-remix
```

- Add the *longranger-[version]* folder to your $PATH so you can run the longranger-remix commands.



## License

The code is licensed under the [GNU Affero General Public License version 3](https://support.10xgenomics.com/docs/license).
