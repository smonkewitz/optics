#pragma once

#include <filesystem>
#include <string_view>

namespace optics {

class InputFile {
    std::string_view data_;

   public:
    explicit InputFile(std::filesystem::path const &path) : InputFile(path.c_str()) {}
    explicit InputFile(char const *path);

    InputFile(InputFile const &) = delete;
    InputFile(InputFile &&) = delete;
    InputFile &operator=(InputFile const &) = delete;
    InputFile &operator=(InputFile &&) = delete;

    ~InputFile();

    std::string_view data() const { return data_; }
};

}  // namespace optics
