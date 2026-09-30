
#pragma once

#include <sstream>
#include <string>
#include <string_view>

enum FileFormat {
    FORMAT_AUTO_DETECT,
    FORMAT_FASTA,
    FORMAT_FASTQ,
    FORMAT_FASTA_QUAL
};

/*enum FastqVersion {
    VERSION_AUTO_DETECT,
    SANGER,
    SOLEXA,
    ILLUMINA_1_3,
    ILLUMINA_1_5,
};*/



struct FastxRecord {
    FileFormat format;
    std::string header = "";  // header line, including @/>, but not newline
    std::string id = "";      // from first char. after @/> up to first whitespace
    std::string sequence = "";
    std::string quality = "";   // only meaningful for FASTQ seqs
    
    std::string &to_string();

private:
    std::string str_representation;
};

class FastxUtils {
public:
    static bool IsValidPair(FastxRecord &a, FastxRecord &b) {
        if (a.id == b.id) return true;
        if (a.id.length() != b.id.length()) return false;
        int i;
        for (i = 0; i < a.id.length(); i++) {
            if (a.id[i] != b.id[i]) break;
        }
        return i == a.id.length() - 1 && a.id[i] == '1' && b.id[i] == '2';
    }
};

class BufferedFastxReader {
private:
    std::stringstream str_stream_;
    std::string str_buffer_ = "";
    FileFormat file_format_;
    char* block_buffer_;
    size_t block_buffer_size_;
    bool strip_space_ = true;
    size_t last_block_size_ = 0;

    // FASTQ batches (LoadBatch): whole lines, each ending in '\n', parsed from batch_pos_ on by
    // NextFastq outside the reader lock. FASTA goes through str_stream_.
    std::string batch_;
    size_t batch_pos_ = 0;
    bool batch_mode_ = false;

    bool m_error = false;

    bool NextBatchLine(std::string_view &line);
    bool NextFastq(FastxRecord &record);
    bool FinishFastq(FastxRecord &record, std::string_view sequence, std::string_view plus, std::string_view quality);

public:
    BufferedFastxReader();
    
    ~BufferedFastxReader();
    
    BufferedFastxReader(const BufferedFastxReader &rhs) = delete;
    
    BufferedFastxReader &operator=(const BufferedFastxReader &rhs) = delete;
    
    bool LoadBatch(std::istream &ifs, size_t record_count);
    
    bool LoadBlock(std::istream &ifs, size_t block_size);
    
    bool NextSequence(FastxRecord &record);

    // static
    bool ReadNextSequence(std::istream &is, FastxRecord &record,
                                 std::string &str_buffer_ptr, FileFormat format = FORMAT_AUTO_DETECT);

    bool Error() {
        return m_error;
    }

    size_t LastBlockSize();

    void auto_detect() { file_format_ = FORMAT_AUTO_DETECT; }
    
    FileFormat file_format() { return file_format_; }
};
