//! 中文：输入读取模块，负责以流式方式解析 FASTA、FASTQ 以及对应的 gzip 压缩文件。
//! English: Streaming input module that parses FASTA, FASTQ, and their gzip-compressed variants.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use flate2::bufread::MultiGzDecoder;

const INPUT_BUFFER_CAPACITY: usize = 1 << 20;

/// 中文：输入记录所属的源格式，用来区分 FASTA 和 FASTQ 解析路径。
/// English: Source format for an input record, used to distinguish FASTA from FASTQ parsing paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceFormat {
    Fasta,
    Fastq,
}

/// 中文：统一的读段记录结构，屏蔽 FASTA/FASTQ 之间的格式差异。
/// English: Unified read record structure that hides the format differences between FASTA and FASTQ.
#[derive(Debug, Clone)]
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
    pub source_format: SourceFormat,
}

/// 中文：面向上层扫描器的流式读取器，内部根据格式分派到具体实现。
/// English: Streaming reader exposed to the scanner, internally dispatching to the format-specific parser.
pub struct RecordReader {
    source_format: SourceFormat,
    inner: ReaderImpl,
    progress: InputProgress,
}

/// 中文：输入进度状态，记录已读取字节数和总字节数，供进度条使用。
/// English: Input progress state tracking bytes read and total bytes for progress reporting.
#[derive(Clone)]
struct InputProgress {
    bytes_read: Arc<AtomicU64>,
    total_bytes: u64,
}

/// 中文：一次进度快照，避免上层直接接触内部原子计数器。
/// English: One immutable progress snapshot so callers do not need to touch the internal atomics directly.
#[derive(Debug, Clone, Copy)]
pub struct ProgressSnapshot {
    pub bytes_read: u64,
    pub total_bytes: u64,
}

/// 中文：内部读取实现枚举，根据检测到的输入类型切换状态机。
/// English: Internal parser enum that switches between format-specific state machines.
enum ReaderImpl {
    Fasta(FastaReader),
    Fastq(FastqReader),
}

/// 中文：打开输入文件并自动检测 FASTA/FASTQ 格式，返回统一的 `RecordReader`。
/// English: Opens the input file, auto-detects FASTA/FASTQ format, and returns a unified `RecordReader`.
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
    /// 中文：返回当前读取器对应的源格式。
    /// English: Returns the source format handled by this reader.
    pub fn source_format(&self) -> SourceFormat {
        self.source_format
    }

    /// 中文：按 chunk 拉取下一批记录；这是并行扫描层的天然批处理边界。
    /// English: Pulls the next chunk of records; this is the natural batch boundary used by the parallel scanner.
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

    /// 中文：获取当前输入读取进度的快照。
    /// English: Returns a snapshot of the current input-reading progress.
    pub fn progress_snapshot(&self) -> ProgressSnapshot {
        self.progress.snapshot()
    }
}

// 中文：打开底层输入流，并在需要时自动套上 gzip 解压和字节计数包装器。
// English: Opens the underlying input stream and wraps it with gzip decoding and byte counting when needed.
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
    // 中文：把内部原子计数转换成稳定快照，同时保证总字节数至少为 1，避免除零问题。
    // English: Converts the internal atomics into a stable snapshot and guarantees total bytes is at least 1 to avoid divide-by-zero issues.
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

/// 中文：为任意底层 reader 增加“读取字节计数”能力的轻量包装器。
/// English: Lightweight wrapper that adds byte-counting behavior to any underlying reader.
struct CountingReader<R> {
    inner: R,
    bytes_read: Arc<AtomicU64>,
}

impl<R> CountingReader<R> {
    /// 中文：创建一个新的计数字节读取器。
    /// English: Creates a new byte-counting reader wrapper.
    fn new(inner: R, bytes_read: Arc<AtomicU64>) -> Self {
        Self { inner, bytes_read }
    }
}

impl<R: Read> Read for CountingReader<R> {
    /// 中文：转发底层 `read`，并把实际读取到的字节数累加到进度计数器。
    /// English: Forwards the underlying `read` call and adds the consumed byte count to the progress counter.
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let bytes = self.inner.read(buf)?;
        self.bytes_read.fetch_add(bytes as u64, Ordering::Relaxed);
        Ok(bytes)
    }
}

impl<R: BufRead> BufRead for CountingReader<R> {
    /// 中文：直接暴露底层缓冲区，不在这里额外复制数据。
    /// English: Exposes the underlying buffer directly without copying data here.
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    /// 中文：消费缓冲区时同步增加已读字节计数。
    /// English: Increments the read-byte counter when the buffer is consumed.
    fn consume(&mut self, amt: usize) {
        self.bytes_read.fetch_add(amt as u64, Ordering::Relaxed);
        self.inner.consume(amt);
    }
}

// 中文：优先根据文件扩展名检测 FASTA/FASTQ；如果扩展名不可靠，再回退到首个非空白字符判断。
// English: Detects FASTA/FASTQ primarily from the file extension and falls back to the first non-whitespace byte when needed.
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

// 中文：根据文件名后缀做快速格式判断，用来减少对文件内容的探测成本。
// English: Performs a fast format guess from the filename suffix to avoid probing file contents when possible.
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

// 中文：读取一行并去掉结尾的换行符，方便 FASTA/FASTQ 状态机直接处理内容。
// English: Reads one line and strips trailing newline characters so the FASTA/FASTQ state machines can process clean content.
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

// 中文：把一行文本序列转成大写字节数组。
// English: Converts one text sequence line into an uppercase byte vector.
fn uppercase_bytes(input: &str) -> Vec<u8> {
    input
        .as_bytes()
        .iter()
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}

// 中文：把去空白并大写后的序列片段追加到目标缓冲区，主要用于 multiline FASTA。
// English: Appends a trimmed, uppercased sequence fragment into the destination buffer, mainly for multiline FASTA parsing.
fn append_uppercase_trimmed(input: &str, sink: &mut Vec<u8>) {
    sink.extend(
        input
            .trim()
            .as_bytes()
            .iter()
            .map(|byte| byte.to_ascii_uppercase()),
    );
}

/// 中文：FASTA 读取器，维护“待处理 header”和当前行缓冲，实现顺序状态解析。
/// English: FASTA reader that maintains a pending header and line buffer to implement sequential state-based parsing.
struct FastaReader {
    reader: Box<dyn BufRead>,
    pending_header: Option<String>,
    line_buffer: String,
    done: bool,
}

impl FastaReader {
    /// 中文：创建 FASTA 读取器并初始化内部状态。
    /// English: Creates a FASTA reader with its internal parser state initialized.
    fn new(reader: Box<dyn BufRead>) -> Self {
        Self {
            reader,
            pending_header: None,
            line_buffer: String::new(),
            done: false,
        }
    }

    /// 中文：读取下一条 FASTA 记录，支持 multiline FASTA，并在遇到下一个 header 时停下。
    /// English: Reads the next FASTA record, supporting multiline FASTA and stopping when the next header is encountered.
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

/// 中文：FASTQ 读取器，按标准四行结构逐条解析记录。
/// English: FASTQ reader that parses records using the standard four-line layout.
struct FastqReader {
    reader: Box<dyn BufRead>,
    line_buffer: String,
}

impl FastqReader {
    /// 中文：创建 FASTQ 读取器。
    /// English: Creates a FASTQ reader.
    fn new(reader: Box<dyn BufRead>) -> Self {
        Self {
            reader,
            line_buffer: String::new(),
        }
    }

    /// 中文：读取下一条 FASTQ 记录，并把质量值从 ASCII Phred+33 转成数值形式。
    /// English: Reads the next FASTQ record and decodes quality scores from ASCII Phred+33 into numeric values.
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
    // 中文：验证 FASTA 读取器能把多行 sequence 正确拼接成一条记录。
    // English: Verifies that the FASTA reader correctly joins multiline sequence content into one record.
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
    // 中文：验证 FASTQ 读取器会把序列转成大写，并把质量值按 Phred+33 解码。
    // English: Verifies that the FASTQ reader uppercases sequence data and decodes qualities as Phred+33 values.
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
    // 中文：验证 gzip 包装输入也能通过统一接口被正确读取。
    // English: Verifies that gzip-compressed input can also be read correctly through the same unified reader interface.
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
