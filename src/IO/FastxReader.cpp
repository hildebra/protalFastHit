

#include "FastxReader.h"
#include "ThreadedGzStream.h"
#include <cctype>
#include <cstring>
#include <err.h>
#include <fcntl.h>
#include <cstdio>
#include <sysexits.h>
#include <iostream>
#include <ctype.h>
;
using std::string;

inline void StripString(string &str) {
    while (isspace(str.back()))
        str.pop_back();
}

inline std::string_view Stripped(std::string_view str) {
    while (!str.empty() && isspace(str.back()))
        str.remove_suffix(1);
    return str;
}

string& FastxRecord::to_string() {
    str_representation.assign(header);
    switch (format) {
        case FORMAT_FASTQ:
            str_representation.append("\n");
            str_representation.append(sequence);
            str_representation.append("\n+\n");
            str_representation.append(quality);
            str_representation.append("\n");
            break;
        default:
            str_representation.append("\n");
            str_representation.append(sequence);
            str_representation.append("\n");
            break;
    }
    return str_representation;
}

BufferedFastxReader::BufferedFastxReader() {
    file_format_ = FORMAT_AUTO_DETECT;
    str_buffer_.reserve(8192);
    block_buffer_ = new char[8192];
    block_buffer_size_ = 8192;
}


BufferedFastxReader::~BufferedFastxReader() {
    delete[] block_buffer_;
}

bool BufferedFastxReader::LoadBlock(std::istream &ifs, size_t block_size) {
    str_stream_.clear();
    str_stream_.str("");
    batch_mode_ = false;
    if (block_buffer_size_ < block_size) {
        delete[] block_buffer_;
        block_buffer_ = new char[block_size];
        block_buffer_size_ = block_size;
    }
    ifs.read(block_buffer_, block_size);
    if (! ifs && ifs.gcount() <= 0)
        return false;
    
    if (file_format_ == FORMAT_AUTO_DETECT) {
        switch (block_buffer_[0]) {
            case '@' : file_format_ = FORMAT_FASTQ; break;
            case '>' : file_format_ = FORMAT_FASTA; break;
            default:
                m_error = true;
                std::cerr << "sequence reader - unrecognized file format " << block_buffer_[0] << std::endl;
                return false;
                errx(EX_DATAERR, "sequence reader - unrecognized file format %c", block_buffer_[0]);
        }
    }

    str_buffer_.assign(block_buffer_, ifs.gcount());
    last_block_size_ = ifs.gcount();
    str_stream_ << str_buffer_;

    if (getline(ifs, str_buffer_)) {
        str_stream_ << str_buffer_ << "\n";
    }
    if (file_format_ == FORMAT_FASTQ) {
        while (getline(ifs, str_buffer_)) {
            str_stream_ << str_buffer_ << "\n";
            if (str_buffer_[0] == '@')
                break;
        }
        int lines_to_read = 0;
        if (getline(ifs, str_buffer_)) {
            str_stream_ << str_buffer_ << "\n";
            lines_to_read = str_buffer_[0] == '@' ? 3 : 2;
            while (lines_to_read-- > 0 && getline(ifs, str_buffer_)) {
                str_stream_ << str_buffer_ << "\n";
            }
        }
    }
    else {
        while (ifs) {
            if (ifs.peek() == '>')
                break;
            if (getline(ifs, str_buffer_))
                str_stream_ << str_buffer_ << "\n";
        }
    }
    return true;
}

bool BufferedFastxReader::LoadBatch(std::istream &ifs, size_t record_count) {
    str_stream_.clear();
    str_stream_.str("");
    batch_.clear();
    batch_pos_ = 0;
    batch_mode_ = false;
    auto valid = false;
    if (file_format_ == FORMAT_AUTO_DETECT) {
        if (!ifs)
            return false;
        switch (ifs.peek()) {
            case '@' :
                file_format_ = FORMAT_FASTQ;
                break;
            case '>' :
                file_format_ = FORMAT_FASTA;
                break;
            case EOF :
                return false;
            default:
                m_error = true;
                std::cerr << "sequence reader - unrecognized file format" << std::endl;
                return false;
                // errx(EX_DATAERR, "sequence reader - unrecognized file format");
        }
        valid = true;
    }

    // FASTQ: the record_count records' 4 lines each, parsed later outside the reader lock
    // (NextFastq). From a ThreadedGzStreambuf, protal's read input, they are cut from its inflated
    // block with memchr and copied at once; from other streams read line by line.
    if (file_format_ == FORMAT_FASTQ) {
        batch_mode_ = true;
        size_t const lines = 4 * record_count;
        size_t taken = 0;
        if (auto* gz = dynamic_cast<protal::ThreadedGzStreambuf*>(ifs.rdbuf())) {
            if (ifs) taken = gz->TakeLines(lines, batch_);
        } else {
            while (taken < lines && getline(ifs, str_buffer_)) {
                batch_.append(str_buffer_);
                batch_.push_back('\n');
                taken++;
            }
        }
        last_block_size_ = batch_.size();
        return valid || taken > 0;
    }

    auto before = ifs.tellg();

    while (record_count > 0 && ifs) {
        str_buffer_.clear();
        if (!getline(ifs, str_buffer_))
            break;
        valid = true;
        if (ifs.peek() == '>')
            record_count--;
        str_stream_ << str_buffer_ << "\n";
    }

    last_block_size_ = ifs.tellg() - before;

    return valid;
}

bool BufferedFastxReader::NextSequence(FastxRecord &seq) {
    if (batch_mode_) return NextFastq(seq);
    return BufferedFastxReader::ReadNextSequence
            (str_stream_, seq, str_buffer_, file_format_);
}

// The next line of batch_, without its '\n'.
bool BufferedFastxReader::NextBatchLine(std::string_view &line) {
    if (batch_pos_ >= batch_.size()) return false;
    char const* begin = batch_.data() + batch_pos_;
    size_t const rest = batch_.size() - batch_pos_;
    auto const* newline = static_cast<char const*>(std::memchr(begin, '\n', rest));
    size_t const length = newline ? static_cast<size_t>(newline - begin) : rest;
    line = std::string_view(begin, length);
    batch_pos_ += length + 1;
    return true;
}

// ReadNextSequence's FASTQ case on the lines of batch_: the same records, ends and errors.
bool BufferedFastxReader::NextFastq(FastxRecord &record) {
    std::string_view header;
    if (!NextBatchLine(header))
        return false;
    header = Stripped(header);
    record.format = FORMAT_FASTQ;
    if (header.empty()) { // Allow empty line to end file
        return false;
    }
    if (header[0] != '@') {
        m_error = true;
        std::cerr << "malformed FASTQ file (exp. '@', saw " << header << ")" << __LINE__ << " in " << __FILE__ << std::endl;
        return false;
    }
    record.header.assign(header);
    if (header.size() <= 1)
        return false;
    auto const first_whitespace_ch = header.find_first_of(" \t\r", 1);
    record.id.assign(header.substr(1, first_whitespace_ch == std::string_view::npos ? std::string_view::npos : first_whitespace_ch - 1));

    std::string_view line;
    if (!NextBatchLine(line))
        return false;
    record.sequence.assign(Stripped(line));
    if (!NextBatchLine(line))  //  + line, discard
        return false;
    if (!NextBatchLine(line))
        return false;
    record.quality.assign(Stripped(line));
    return true;
}

bool BufferedFastxReader::ReadNextSequence(std::istream &is, FastxRecord &record,
                                           std::string &str_buffer, FileFormat file_format) {
    str_buffer.clear();

    if (!getline(is, str_buffer)) {
        return false;
    }

    StripString(str_buffer);
    if (file_format == FORMAT_AUTO_DETECT) {
        switch (str_buffer[0]) {
            case '@' :
                file_format = FORMAT_FASTQ;
                break;
            case '>' :
                if (!isdigit(is.peek())) {
                    file_format = FORMAT_FASTA;
                } else {
                    file_format = FORMAT_FASTA_QUAL;
                }
                break;
            default:
                m_error = true;
                std::cerr << "sequence reader - unrecognized file format" << std::endl;
                return false;
                // errx(EX_DATAERR, "sequence reader - unrecognized file format");
        }
    }
    record.format = file_format;
    if (record.format == FORMAT_FASTQ) {
        if (str_buffer.empty()) { // Allow empty line to end file
            return false;
        }
        if (str_buffer[0] != '@') {
            m_error = true;
            std::cerr << "malformed FASTQ file (exp. '@', saw " << str_buffer.c_str() << ")" << __LINE__ << " in " << __FILE__ << std::endl;
            return false;
            // errx(EX_DATAERR, "malformed FASTQ file (exp. '@', saw \"%s\"), aborting",
            //      str_buffer.c_str());
        }
    } else if (record.format == FORMAT_FASTA || record.format == FORMAT_FASTA_QUAL) {
        if (str_buffer[0] != '>') {
            m_error = true;
            std::cerr << "malformed FASTQ file (exp. '>', saw " << str_buffer.c_str() << ")" << __LINE__ << " in " << __FILE__ << std::endl;
            return false;
            // errx(EX_DATAERR, "malformed FASTA file (exp. '>', saw \"%s\"), aborting",
            //      str_buffer.c_str());
        }
    } else {
        m_error = true;
        std::cerr << "illegal sequence format encountered in parsing" << std::endl;
        return false;
        // errx(EX_SOFTWARE, "illegal sequence format encountered in parsing");
    }
    record.header.assign(str_buffer);
    auto first_whitespace_ch = str_buffer.find_first_of(" \t\r", 1);
    auto substr_len = first_whitespace_ch;
    if (substr_len != std::string::npos)
        substr_len--;
    if (str_buffer.size() > 1)
        record.id.assign(str_buffer, 1, substr_len);
    else {
        return false;
    }

    if (record.format == FORMAT_FASTQ) {
        if (!getline(is, str_buffer))
            return false;
        StripString(str_buffer);
        record.sequence.assign(str_buffer);
        if (!getline(is, str_buffer))  //  + line, discard
            return false;
        if (!getline(is, str_buffer))
            return false;
        StripString(str_buffer);
        record.quality.assign(str_buffer);

    } else if (record.format == FORMAT_FASTA) {
        record.quality.assign("");
        record.sequence.assign("");
        while (is && is.peek() != '>') {
            if (!getline(is, str_buffer))
                return !record.sequence.empty();
            StripString(str_buffer);
            record.sequence.append(str_buffer);
        }
    } else if (record.format == FORMAT_FASTA_QUAL) {
        record.quality.assign("");
        record.sequence.assign("");
        while (is && is.peek() != '>') {
            if (!getline(is, str_buffer))
                return !record.sequence.empty();
            StripString(str_buffer);
            record.sequence.append(str_buffer);
            record.sequence.append(" ");
        }
        StripString(record.sequence);
    }
    return true;
}

size_t BufferedFastxReader::LastBlockSize() {
    return last_block_size_;
}

