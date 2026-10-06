// make_journal.cpp - writes the journal that ReplaceTail would write for OLD_DB (a copy of a database before an
// in-place --add_model) as DB.journal, from its first model member on: with DB cut short, it is a stopped
// --add_model.
//
//   make_journal OLD_DB DB
#include "Utilities/Database.h"
#include "ReadType.h"

#include <iostream>

using namespace protal;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: make_journal OLD_DB DB" << std::endl;
        return 1;
    }
    std::string error;
    auto const bundle = db::Bundle::Open(argv[1], error);
    if (!bundle) {
        std::cerr << argv[1] << ": " << error << std::endl;
        return 1;
    }
    auto const models = AllModelFiles();
    uint64_t tail = 0;
    for (auto const& member : bundle->Members()) {
        if (std::find(models.begin(), models.end(), member.name) != models.end()) {
            tail = member.frames.frames.front().compressed_offset;
            break;
        }
    }
    std::ifstream is(argv[1], std::ios::binary);
    std::string const bytes{std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>()};
    uint64_t const directory = bundle->DirectoryFrame().compressed_size;
    db::detail::Journal journal;
    journal.tail_offset = tail;
    journal.check_offset = std::max<uint64_t>(directory, tail >= 4096 ? tail - 4096 : 0);
    journal.directory = bytes.substr(0, directory);
    journal.check = bytes.substr(journal.check_offset, tail - journal.check_offset);
    journal.tail = bytes.substr(tail);
    std::ofstream(std::string(argv[2]) + db::kJournalExtension, std::ios::binary) << journal.Serialize();
    std::cout << "journal: tail at " << tail << ", " << journal.tail.size() << " bytes" << std::endl;
    return 0;
}
