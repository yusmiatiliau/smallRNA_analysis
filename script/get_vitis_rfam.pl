#!/usr/bin/env perl

use strict;
use warnings;

use DBI;
use File::Path qw(make_path);

# Rfam public MySQL database
my $host     = 'mysql-rfam-public.ebi.ac.uk';
my $port     = 4497;
my $database = 'Rfam';
my $user     = 'rfamro';
my $password = '';

# NCBI taxonomy ID for Vitis vinifera
my $taxid = 29760;

my $output_dir = 'vitis_rfam_output';

make_path($output_dir);
make_path("$output_dir/regions_by_family");

my $dsn = join(
    ';',
    "DBI:mysql:database=$database",
    "host=$host",
    "port=$port"
);

print STDERR "Connecting to the Rfam database...\n";

my $dbh = DBI->connect(
    $dsn,
    $user,
    $password,
    {
        RaiseError => 1,
        PrintError => 0,
        AutoCommit => 1,
    }
) or die "Could not connect to Rfam: $DBI::errstr\n";

my $sql = q{
    SELECT DISTINCT
           CASE
               WHEN f.type LIKE '%tRNA%' THEN 'tRNA'
               WHEN f.type LIKE '%rRNA%' THEN 'rRNA'
           END AS rna_class,
           fr.rfam_acc,
           CONCAT(
               fr.rfamseq_acc,
               '/',
               fr.seq_start,
               '-',
               fr.seq_end
           ) AS sequence_region
    FROM full_region AS fr
    JOIN rfamseq AS rs
      ON rs.rfamseq_acc = fr.rfamseq_acc
    JOIN family AS f
      ON f.rfam_acc = fr.rfam_acc
    WHERE rs.ncbi_id = ?
      AND (
            f.type LIKE '%tRNA%'
            OR f.type LIKE '%rRNA%'
          )
      AND fr.is_significant = 1
      AND fr.type = 'full'
    ORDER BY rna_class, fr.rfam_acc, sequence_region
};

my $sth = $dbh->prepare($sql);
$sth->execute($taxid);

open my $table_fh, '>',
    "$output_dir/Vitis_vinifera_tRNA_rRNA_by_family.tsv"
    or die "Cannot write output table: $!\n";

print {$table_fh}
    "rna_class\trfam_accession\tsequence_region\n";

my %class_regions;
my %family_regions;

while (my ($rna_class, $rfam_acc, $region) =
       $sth->fetchrow_array()) {

    next unless defined $rna_class;
    next unless defined $rfam_acc;
    next unless defined $region;

    print {$table_fh}
        join("\t", $rna_class, $rfam_acc, $region), "\n";

    $class_regions{$rna_class}{$region} = 1;
    $family_regions{$rna_class}{$rfam_acc}{$region} = 1;
}

close $table_fh;

$sth->finish();
$dbh->disconnect();

# Write combined coordinate files for tRNA and rRNA.
for my $rna_class (qw(tRNA rRNA)) {

    my $region_file =
        "$output_dir/Vitis_vinifera_${rna_class}_regions.txt";

    open my $fh, '>', $region_file
        or die "Cannot write $region_file: $!\n";

    my @regions = sort keys %{
        $class_regions{$rna_class} // {}
    };

    for my $region (@regions) {
        print {$fh} "$region\n";
    }

    close $fh;

    print STDERR
        "$rna_class regions: ",
        scalar(@regions),
        "\n";
}

# Write one coordinate file for each Rfam family.
for my $rna_class (sort keys %family_regions) {

    for my $rfam_acc (
        sort keys %{ $family_regions{$rna_class} }
    ) {

        my $region_file =
            "$output_dir/regions_by_family/" .
            "${rna_class}_${rfam_acc}.txt";

        open my $fh, '>', $region_file
            or die "Cannot write $region_file: $!\n";

        for my $region (
            sort keys %{
                $family_regions{$rna_class}{$rfam_acc}
            }
        ) {
            print {$fh} "$region\n";
        }

        close $fh;
    }
}

print STDERR "Results written to $output_dir\n";