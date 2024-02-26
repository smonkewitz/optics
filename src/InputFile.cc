#include "InputFile.h"

#include <absl/cleanup/cleanup.h>
#include <fcntl.h>
#include <fmt/core.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cerrno>
#include <stdexcept>

namespace optics {

InputFile::InputFile(char const *path) {
    int fd = ::open(path, O_RDONLY | O_CLOEXEC);
    if (fd == -1) {
        throw std::runtime_error(
            fmt::format("failed to open {}: errno={}", path, errno));
    }

    absl::Cleanup const closer = [fd] { ::close(fd); };
    struct ::stat buf;
    std::memset(&buf, 0, sizeof(struct ::stat));
    if (fstat(fd, &buf) == -1) {
        throw std::runtime_error(
            fmt::format("failed to fstat {}: errno={}", path, errno));
    }
    if (buf.st_size <= 0 || buf.st_size > std::numeric_limits<ssize_t>::max()) {
        throw std::runtime_error("File is empty, or too large to map into memory");
    }
    auto const size = static_cast<size_t>(buf.st_size);
    int flags = MAP_PRIVATE | MAP_NORESERVE;
#if defined(MAP_HUGETLB) && defined(MAP_HUGE_2MB) && defined(MAP_POPULATE)
    flags |= MAP_HUGETLB | MAP_HUGE_2MB | MAP_POPULATE;
#endif
    void *data = ::mmap(nullptr, size, PROT_READ, flags, fd, 0);
    if (data == MAP_FAILED) {
        throw std::runtime_error(
            fmt::format("failed to mmap contents of {}: errno={}", path, errno));
    }
    if (::madvise(data, size, MADV_WILLNEED) == -1) {
        throw std::runtime_error(
            fmt::format("madvise on contents of {} failed: errno={}", path, errno));
    }
    data_ = std::string_view{static_cast<const char *>(data), size};
}

InputFile::~InputFile() { ::munmap(const_cast<char *>(data_.data()), data_.size()); }

}  // namespace optics
