#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKDIR="${BENCH_WORKDIR:-$ROOT_DIR/benchmark_work}"
READ_COUNT="${BENCH_READS:-20000}"
MOTIF_COUNT="${BENCH_MOTIFS:-128}"
READ_LEN="${BENCH_READ_LEN:-160}"
SEED="${BENCH_SEED:-1337}"
INSERT_POS="${BENCH_INSERT_POS:-20}"
PACKAGE_NAME="motifscan-benchmark"
PACKAGE_DIR="$WORKDIR/package"
INSTALL_DIR="$WORKDIR/install"
DATA_DIR="$WORKDIR/data"
RESULT_DIR="$WORKDIR/results"
ARCHIVE_NAME="$PACKAGE_NAME-$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m).tar.gz"

mkdir -p "$WORKDIR"
WORKDIR_ABS="$(cd "$WORKDIR" && pwd)"
case "$WORKDIR_ABS" in
  "$ROOT_DIR"/*) ;;
  *)
    echo "refusing to use a workdir outside the repository: $WORKDIR_ABS" >&2
    exit 1
    ;;
esac

if [[ "$WORKDIR_ABS" == "$ROOT_DIR" ]]; then
  echo "refusing to use the repository root as the benchmark workdir" >&2
  exit 1
fi

rm -rf "$WORKDIR_ABS"
WORKDIR="$WORKDIR_ABS"
mkdir -p "$PACKAGE_DIR" "$INSTALL_DIR" "$DATA_DIR" "$RESULT_DIR"

time_step() {
  local label="$1"
  shift
  local start_ns end_ns elapsed_ns elapsed_s
  start_ns="$(date +%s%N)"
  "$@"
  end_ns="$(date +%s%N)"
  elapsed_ns=$((end_ns - start_ns))
  elapsed_s="$(awk "BEGIN { printf \"%.3f\", $elapsed_ns / 1000000000 }")"
  printf '%s: %s s\n' "$label" "$elapsed_s"
}

compare_files() {
  local expected="$1"
  local actual="$2"
  local label="$3"
  if cmp -s "$expected" "$actual"; then
    printf '%s: ok\n' "$label"
  else
    printf '%s: mismatch\n' "$label" >&2
    diff -u "$expected" "$actual" || true
    exit 1
  fi
}

generate_data() {
  if (( READ_LEN <= INSERT_POS + 16 )); then
    echo "READ_LEN must be greater than INSERT_POS + 16" >&2
    exit 1
  fi

  perl - "$DATA_DIR" "$READ_COUNT" "$MOTIF_COUNT" "$READ_LEN" "$SEED" "$INSERT_POS" <<'PERL'
use strict;
use warnings;

my ($dir, $read_count, $motif_count, $read_len, $seed, $insert_pos) = @ARGV;
srand($seed);

sub rc {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGT/TGCA/;
    return $seq;
}

sub rand_base {
    my @bases = qw(A C G T);
    return $bases[int(rand(@bases))];
}

sub random_motif {
    my $len = 6 + int(rand(9));
    my $seq = 'A';
    for (1 .. $len - 1) {
        $seq .= rand_base();
    }
    return $seq;
}

my %seen;
my @motifs;
while (@motifs < $motif_count) {
    my $candidate = random_motif();
    my $candidate_rc = rc($candidate);
    next if $candidate eq $candidate_rc;
    next if exists $seen{$candidate};
    next if exists $seen{$candidate_rc};
    $seen{$candidate} = 1;
    $seen{$candidate_rc} = 1;
    push @motifs, $candidate;
}

open my $motif_fh, '>', "$dir/motifs.csv" or die $!;
print {$motif_fh} "name,sequence\n";
for my $i (0 .. $#motifs) {
    print {$motif_fh} 'motif', ($i + 1), ',', $motifs[$i], "\n";
}
close $motif_fh;

open my $fq_fh, '>', "$dir/reads.fastq" or die $!;
open my $fa_fh, '>', "$dir/reads.fa" or die $!;
open my $summary_fh, '>', "$dir/expected_summary.csv" or die $!;
open my $hits_fh, '>', "$dir/expected_hits.csv" or die $!;

print {$summary_fh} "motif,sequence,length,reads_with_hit,total_hits,forward_hits,revcomp_hits\n";
print {$hits_fh} "read_id,motif,strand,position,matched_sequence\n";

my @stats;
for my $motif (@motifs) {
    push @stats, {
        reads_with_hit => 0,
        total_hits => 0,
        forward_hits => 0,
        revcomp_hits => 0,
    };
}

for my $i (1 .. $read_count) {
    my $motif_index = ($i - 1) % @motifs;
    my $motif_name = 'motif' . ($motif_index + 1);
    my $motif = $motifs[$motif_index];
    my $strand = ($i % 2 == 0) ? 'forward' : 'revcomp';
    my $inserted = $strand eq 'forward' ? $motif : rc($motif);

    my $background = ('CGT' x (($read_len / 3) + 3));
    my $seq = substr($background, 0, $read_len);
    substr($seq, $insert_pos, length($inserted)) = $inserted;

    my $read_id = "r$i";
    print {$fq_fh} "\@$read_id\n$seq\n+\n", ('I' x length($seq)), "\n";
    print {$fa_fh} ">$read_id\n$seq\n";

    $stats[$motif_index]{reads_with_hit}++;
    $stats[$motif_index]{total_hits}++;
    if ($strand eq 'forward') {
        $stats[$motif_index]{forward_hits}++;
    } else {
        $stats[$motif_index]{revcomp_hits}++;
    }

    print {$hits_fh} join(',', $read_id, $motif_name, $strand, $insert_pos, $inserted), "\n";
}

for my $i (0 .. $#motifs) {
    my $motif_name = 'motif' . ($i + 1);
    my $stats = $stats[$i];
    print {$summary_fh} join(
        ',',
        $motif_name,
        $motifs[$i],
        length($motifs[$i]),
        $stats->{reads_with_hit},
        $stats->{total_hits},
        $stats->{forward_hits},
        $stats->{revcomp_hits}
    ), "\n";
}

close $fq_fh;
close $fa_fh;
close $summary_fh;
close $hits_fh;
PERL
}

build_release() {
  cargo build --release
  cp "$ROOT_DIR/target/release/motifscan" "$PACKAGE_DIR/motifscan"
  tar -czf "$WORKDIR/$ARCHIVE_NAME" -C "$PACKAGE_DIR" motifscan
  tar -xzf "$WORKDIR/$ARCHIVE_NAME" -C "$INSTALL_DIR"
}

run_benchmark() {
  local bin="$INSTALL_DIR/motifscan"
  local fastq_summary="$RESULT_DIR/count.fastq.csv"
  local fasta_summary="$RESULT_DIR/count.fasta.csv"
  local read_hits="$RESULT_DIR/read_hits.csv"

  time_step "benchmark run (FASTQ + read hits)" "$bin" count \
    -i "$DATA_DIR/reads.fastq" \
    --motifs "$DATA_DIR/motifs.csv" \
    --revcomp \
    --report-read-hits "$read_hits" \
    -o "$fastq_summary"

  compare_files "$DATA_DIR/expected_summary.csv" "$fastq_summary" "FASTQ summary"
  compare_files "$DATA_DIR/expected_hits.csv" "$read_hits" "read-hit output"

  time_step "benchmark run (FASTA summary)" "$bin" count \
    -i "$DATA_DIR/reads.fa" \
    --motifs "$DATA_DIR/motifs.csv" \
    --revcomp \
    -o "$fasta_summary"

  compare_files "$DATA_DIR/expected_summary.csv" "$fasta_summary" "FASTA summary"
}

printf 'Workdir: %s\n' "$WORKDIR"
printf 'Reads: %s\n' "$READ_COUNT"
printf 'Motifs: %s\n' "$MOTIF_COUNT"
printf 'Read length: %s\n' "$READ_LEN"
printf 'Seed: %s\n' "$SEED"

time_step "release build and packaging" build_release
time_step "data generation" generate_data
time_step "benchmark and verification" run_benchmark

printf 'Archive: %s\n' "$WORKDIR/$ARCHIVE_NAME"
printf 'Results: %s\n' "$RESULT_DIR"
