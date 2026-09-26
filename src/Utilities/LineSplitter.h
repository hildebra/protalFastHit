//
// Created by fritsche on 18/10/22.
//

#pragma once

#include <vector>
#include <string>
#include <stdexcept>

namespace protal {
    class LineSplitter {
        std::vector<std::string> m_tokens;
        std::string m_delimiter = "\t";

    public:
        // Splits line at every occurrence of delimiter. Empty fields are kept ("a\t\tb" gives "a", "", "b"
        // and "a\t" gives "a", ""), so columns keep their positions; an empty line gives no tokens.
        static void Split(std::string const& line, std::string const& delimiter, std::vector<std::string> &tokens) {
            tokens.clear();
            if (line.empty()) return;
            if (delimiter.empty()) {
                tokens.emplace_back(line);
                return;
            }
            size_t start = 0;
            for (size_t pos; (pos = line.find(delimiter, start)) != std::string::npos; start = pos + delimiter.size()) {
                tokens.emplace_back(line, start, pos - start);
            }
            tokens.emplace_back(line, start);
        }

        std::vector<std::string>& Tokens() {
            return m_tokens;
        }

        void Split(std::string const& line) {
            Split(line, m_delimiter, m_tokens);
        }

        std::string operator[] (size_t index) {
            if (index >= m_tokens.size()) {
                throw std::out_of_range("Line splitter: trying to access index " + std::to_string(index) + " of split line. Number of tokens is " + std::to_string(m_tokens.size()));
            }
            return m_tokens[index];
        }
    };
}
