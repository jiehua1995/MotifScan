use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use flate2::bufread::MultiGzDecoder;

const INPUT_BUFFER_CAPACITY: usize = 1 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone)]
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
    pub source_format: SourceFormat,
}

pub struct RecordReader {
    source_format: SourceFormat,
    inner: ReaderImpl,
    progress: InputProgress,
}

#[derive(Clone)]
struct InputProgress {
    bytes_read: Arc<AtomicU64>,
    total_bytes: u64,
}

#[derive(Debug, Clone, Copy)]
pub struct ProgressSnapshot {
    pub bytes_read: u64,
    pub total_bytes: u64,
}

enum ReaderImpl {
    Fasta(FastaReader),
    Fastq(FastqReader),
}

pub fn open_record_reader(path: &Path) -> Result<RecordReader> {
    let (mut reader, progress) = open_bufread(path)?;
    let source_format = detect_source_format(path, reader.as_mut())?;
    let inner = match source_format {
        SourceFormat::Fasta => ReaderImpl::Fasta(FastaReader::new(reader)),
        SourceFormat::Fastq => ReaderImpl::Fastq(FastqReader::new(reader)),
    };

    Ok(RecordReader {
        source_format,
        inner,
        progress,
    })
}

impl RecordReader {
    pub fn source_format(&self) -> SourceFormat {
        self.source_format
    }

    pub fn next_chunk(&mut self, chunk_size: usize) -> Result<Vec<Record>> {
        let mut records = Vec::with_capacity(chunk_size);
        for _ in 0..chunk_size {
            let next = match &mut self.inner {
                ReaderImpl::Fasta(reader) => reader.next_record()?,
                ReaderImpl::Fastq(reader) => reader.next_record()?,
            };
            match next {
                Some(record) => records.push(record),
                None => break,
            }
        }
        Ok(records)
    }

    pub fn progress_snapshot(&self) -> ProgressSnapshot {
        self.progress.snapshot()
    }
}

fn open_bufread(path: &Path) -> Result<(Box<dyn BufRead>, InputProgress)> {
    let file =
        File::open(path).with_context(|| format!("failed to open input {}", path.display()))?;
    let total_bytes = file
        .metadata()
        .with_context(|| format!("failed to stat input {}", path.display()))?
        .len();
    let progress = InputProgress {
        bytes_read: Arc::new(AtomicU64::new(0)),
        total_bytes,
    };
    let buffered_file = CountingReader::new(
        BufReader::with_capacity(INPUT_BUFFER_CAPACITY, file),
        progress.bytes_read.clone(),
    );
    if path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
    {
        let decoder = MultiGzDecoder::new(buffered_file);
        Ok((
            Box::new(BufReader::with_capacity(INPUT_BUFFER_CAPACITY, decoder)),
            progress,
        ))
    } else {
        Ok((Box::new(buffered_file), progress))
    }
}

impl InputProgress {
    fn snapshot(&self) -> ProgressSnapshot {
        ProgressSnapshot {
            bytes_read: self
                .bytes_read
                .load(Ordering::Relaxed)
                .min(self.total_bytes),
            total_bytes: self.total_bytes.max(1),
        }
    }
}

struct CountingReader<R> {
    inner: R,
    bytes_read: Arc<AtomicU64>,
}

impl<R> CountingReader<R> {
    fn new(inner: R, bytes_read: Arc<AtomicU64>) -> Self {
        Self { inner, bytes_read }
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let bytes = self.inner.read(buf)?;
        self.bytes_read.fetch_add(bytes as u64, Ordering::Relaxed);
        Ok(bytes)
    }
}

impl<R: BufRead> BufRead for CountingReader<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.bytes_read.fetch_add(amt as u64, Ordering::Relaxed);
        self.inner.consume(amt);
    }
}

fn detect_source_format(path: &Path, reader: &mut dyn BufRead) -> Result<SourceFormat> {
    if let Some(format) = format_from_extension(path) {
        return Ok(format);
    }

    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            bail!("input appears to be empty")
        }

        if let Some((_, byte)) = buffer
            .iter()
            .enumerate()
            .find(|(_, byte)| !byte.is_ascii_whitespace())
        {
            if *byte == b'>' {
                return Ok(SourceFormat::Fasta);
            }
            if *byte == b'@' {
                return Ok(SourceFormat::Fastq);
            }
            bail!("could not detect FASTA/FASTQ format from input")
        }

        let consume_len = buffer.len();
        reader.consume(consume_len);
    }
}

fn format_from_extension(path: &Path) -> Option<SourceFormat> {
    let file_name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    if file_name.ends_with(".fa")
        || file_name.ends_with(".fasta")
        || file_name.ends_with(".fa.gz")
        || file_name.ends_with(".fasta.gz")
    {
        return Some(SourceFormat::Fasta);
    }
    if file_name.ends_with(".fq")
        || file_name.ends_with(".fastq")
        || file_name.ends_with(".fq.gz")
        || file_name.ends_with(".fastq.gz")
    {
        return Some(SourceFormat::Fastq);
    }
    None
}

fn read_line_trimmed(reader: &mut dyn BufRead, buffer: &mut String) -> Result<usize> {
    buffer.clear();
    let bytes = reader.read_line(buffer)?;
    if bytes == 0 {
        return Ok(0);
    }
    while matches!(buffer.chars().last(), Some('\n' | '\r')) {
        buffer.pop();
    }
    Ok(bytes)
}

fn uppercase_bytes(input: &str) -> Vec<u8> {
    input
        .as_bytes()
        .iter()
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}

fn append_uppercase_trimmed(input: &str, sink: &mut Vec<u8>) {
    sink.extend(
        input
            .trim()
            .as_bytes()
            .iter()
            .map(|byte| byte.to_ascii_uppercase()),
    );
}

struct FastaReader {
    reader: Box<dyn BufRead>,
    pending_header: Option<String>,
    line_buffer: String,
    done: bool,
}

impl FastaReader {
    fn new(reader: Box<dyn BufRead>) -> Self {
        Self {
            reader,
            pending_header: None,
            line_buffer: String::new(),
            done: false,
        }
    }

    fn next_record(&mut self) -> Result<Option<Record>> {
        if self.done {
            return Ok(None);
        }

        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                let bytes = read_line_trimmed(self.reader.as_mut(), &mut self.line_buffer)?;
                if bytes == 0 {
                    self.done = true;
                    return Ok(None);
                }
                if self.line_buffer.is_empty() {
                    continue;
                }
                if let Some(rest) = self.line_buffer.strip_prefix('>') {
                    break rest.trim().to_string();
                }
                bail!("invalid FASTA: expected header line starting with '>'")
            }
        };

        let mut sequence = Vec::new();
        loop {
            let bytes = read_line_trimmed(self.reader.as_mut(), &mut self.line_buffer)?;
            if bytes == 0 {
                self.done = true;
                break;
            }
            if let Some(rest) = self.line_buffer.strip_prefix('>') {
                self.pending_header = Some(rest.trim().to_string());
                break;
            }
            if !self.line_buffer.is_empty() {
                append_uppercase_trimmed(&self.line_buffer, &mut sequence);
            }
        }

        if sequence.is_empty() {
            bail!("FASTA record '{header}' has an empty sequence")
        }

        Ok(Some(Record {
            id: header,
            seq: sequence,
            qual: None,
            source_format: SourceFormat::Fasta,
        }))
    }
}

struct FastqReader {
    reader: Box<dyn BufRead>,
    line_buffer: String,
}

impl FastqReader {
    fn new(reader: Box<dyn BufRead>) -> Self {
        Self {
            reader,
            line_buffer: String::new(),
        }
    }

    fn next_record(&mut self) -> Result<Option<Record>> {
        let bytes = loop {
            let bytes = read_line_trimmed(self.reader.as_mut(), &mut self.line_buffer)?;
            if bytes == 0 {
                return Ok(None);
            }
            if !self.line_buffer.is_empty() {
                break bytes;
            }
        };

        if bytes == 0 {
            return Ok(None);
        }
        let header = self
            .line_buffer
            .strip_prefix('@')
            .ok_or_else(|| anyhow!("invalid FASTQ: expected header line starting with '@'"))?
            .trim()
            .to_string();

        read_line_trimmed(self.reader.as_mut(), &mut self.line_buffer)?;
        if self.line_buffer.is_empty() {
            bail!("FASTQ record '{header}' has an empty sequence")
        }
        let seq = uppercase_bytes(self.line_buffer.trim());

        read_line_trimmed(self.reader.as_mut(), &mut self.line_buffer)?;
        if !self.line_buffer.starts_with('+') {
            bail!("invalid FASTQ: expected '+' separator for record '{header}'")
        }

        read_line_trimmed(self.reader.as_mut(), &mut self.line_buffer)?;
        if self.line_buffer.len() != seq.len() {
            bail!(
                "FASTQ quality length mismatch for record '{header}': expected {}, got {}",
                seq.len(),
                self.line_buffer.len()
            )
        }
        let qual = self
            .line_buffer
            .as_bytes()
            .iter()
            .map(|byte| {
                if *byte < 33 {
                    Err(anyhow!("invalid FASTQ quality byte in record '{header}'"))
                } else {
                    Ok(byte - 33)
                }
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(Record {
            id: header,
            seq,
            qual: Some(qual),
            source_format: SourceFormat::Fastq,
        }))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::write::GzEncoder;
    use flate2::Compression;
    use tempfile::tempdir;

    use super::{open_record_reader, SourceFormat};

    #[test]
    fn parses_multiline_fasta() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("input.fa");
        std::fs::write(&path, ">r1\nACG\nTT\n>r2\nGGAA\n").unwrap();

        let mut reader = open_record_reader(&path).unwrap();
        assert_eq!(reader.source_format(), SourceFormat::Fasta);

        let chunk = reader.next_chunk(10).unwrap();
        assert_eq!(chunk.len(), 2);
        assert_eq!(chunk[0].id, "r1");
        assert_eq!(chunk[0].seq, b"ACGTT");
        assert_eq!(chunk[1].seq, b"GGAA");
    }

    #[test]
    fn parses_fastq_and_decodes_phred33() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("input.fastq");
        std::fs::write(&path, "@r1\nacgt\n+\nIIII\n").unwrap();

        let mut reader = open_record_reader(&path).unwrap();
        assert_eq!(reader.source_format(), SourceFormat::Fastq);
        let chunk = reader.next_chunk(10).unwrap();
        assert_eq!(chunk.len(), 1);
        assert_eq!(chunk[0].seq, b"ACGT");
        assert_eq!(chunk[0].qual.as_ref().unwrap(), &vec![40, 40, 40, 40]);
    }

    #[test]
    fn parses_gzip_smoke() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("input.fastq.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();

        let mut reader = open_record_reader(&path).unwrap();
        let chunk = reader.next_chunk(10).unwrap();
        assert_eq!(chunk.len(), 1);
        assert_eq!(chunk[0].id, "r1");
    }
}
