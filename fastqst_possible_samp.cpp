#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <sstream>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <optional>
#include <queue>

#include <zlib.h>


template<class T>
class BlockingQueue
{
public:
  explicit BlockingQueue(size_t max_size) : max_size_(max_size) {}

  void enqueue(T obj)
  {
    {
      std::unique_lock<std::mutex> lock(mu_);
      not_full_.wait(lock, [this] { return q_.size() < max_size_; });
      q_.push(std::move(obj));
    }
    not_empty_.notify_one();
  }

  std::optional<T> dequeue()
  {
    std::unique_lock<std::mutex> lock(mu_);
    not_empty_.wait(lock, [this] { return !q_.empty() || shutdown_; });
    if (shutdown_)
      return std::nullopt;
    T obj = std::move(q_.front());
    q_.pop();
    lock.unlock();
    not_full_.notify_one();
    return std::move(obj);
  }

  void shutdown()
  {
    {
      std::unique_lock<std::mutex> lock(mu_);
      shutdown_ = true;
    }
    not_empty_.notify_one();
  }

private:
  const size_t max_size_;

  std::queue<T> q_;
  bool shutdown_ = false;
  std::mutex mu_;
  std::condition_variable not_empty_;
  std::condition_variable not_full_;
};


class GzipReader
{
private:
  gzFile file = nullptr;
  char buffer[4096];
  std::string leftover;

public:
  explicit GzipReader(std::string path)
  {
    file = gzopen(path.c_str(), "rb");
    if (!file)
      throw std::runtime_error("Could not open gzip file: " + path);
  }

  ~GzipReader()
  {
    if (file)
      gzclose(file);
  }

  // includes trailing newline
  bool getline(std::string& line)
  {
    return _getline(line, true);
  }
  bool _getline(std::string& line, bool actually_get)
  {
    line.clear();

    // First use any leftover data from previous read
    if (!leftover.empty())
    {
      size_t pos = leftover.find('\n');
      if (pos != std::string::npos)
      {
        if (actually_get)
          line = leftover.substr(0, pos + 1);
        leftover = leftover.substr(pos + 1);
        return true;
      }
      if (actually_get)
        line = leftover;
      leftover.clear();
    }

    while (true)
    {
      int bytes_read = gzread(file, buffer, sizeof(buffer));
      if (bytes_read <= 0)
        return !line.empty(); // file ends in broken partial line; allow it

      for (int i = 0; i < bytes_read; i++)
      {
        if (buffer[i] == '\n')
        {
          if (actually_get)
            line.append(buffer, i+1);
          leftover.assign(buffer + i + 1, bytes_read - (i + 1));
          return true;
        }
      }
      if (actually_get)
        line.append(buffer, bytes_read);
    }
  }
  std::string junk;
  void skip4lines()
  {
    _getline(junk, false);
    _getline(junk, false);
    _getline(junk, false);
    _getline(junk, false);
  }
};

struct FastqHalfRead
{
  std::string lines[4];
  explicit FastqHalfRead(GzipReader& gzr)
  {
    for (int i = 0; i < 4; i++)
      gzr.getline(lines[i]);
  }
};

class SubsampleWriter
{
public:
  std::vector<uint32_t> selected_indices;
  uint32_t current_index1 = 0; // index into selected_indices
  uint32_t current_index2 = 0; // index into selected_indices
  std::string r1_path;
  std::string r2_path;
  BlockingQueue<FastqHalfRead> q1;
  BlockingQueue<FastqHalfRead> q2;

  gzFile r1_out = nullptr;
  gzFile r2_out = nullptr;

  SubsampleWriter(std::string path1, std::string path2)
  : r1_path(path1), r2_path(path2), q1(1024), q2(1024)
  {
    r1_out = gzopen(r1_path.c_str(), "wb");
    if (!r1_out)
      throw std::runtime_error("Could not open gzip file for writing: " + r1_path);
    r2_out = gzopen(r2_path.c_str(), "wb");
    if (!r2_out)
      throw std::runtime_error("Could not open gzip file for writing: " + r2_path);
  }
  ~SubsampleWriter()
  {
    if (r1_out)
      gzclose(r1_out);
    if (r2_out)
      gzclose(r2_out);
  }

  bool wantThisRead1(uint32_t the_read)
  {
    bool ret = current_index1 < selected_indices.size() &&
           selected_indices[current_index1] == the_read;
    if (ret)
      current_index1++;
    return ret;
  }
  bool wantThisRead2(uint32_t the_read)
  {
    bool ret = current_index2 < selected_indices.size() &&
           selected_indices[current_index2] == the_read;
    if (ret)
      current_index2++;
    return ret;
  }

  void outputThread1()
  {
    while (true)
    {
      auto maybe_read = q1.dequeue();
      if (!maybe_read.has_value())
        break;
      FastqHalfRead const& our_read = maybe_read.value();
      for (int i = 0; i < 4; i++)
        if (gzwrite(r1_out, our_read.lines[i].c_str(), our_read.lines[i].length()) < 1)
        {
          std::cerr << "Failed to write to gzip file " << r1_path << std::endl;
          return;
        }
    }
  }
  void outputThread2()
  {
    while (true)
    {
      auto maybe_read = q2.dequeue();
      if (!maybe_read.has_value())
        break;
      FastqHalfRead const& our_read = maybe_read.value();
      for (int i = 0; i < 4; i++)
        if (gzwrite(r2_out, our_read.lines[i].c_str(), our_read.lines[i].length()) < 1)
        {
          std::cerr << "Failed to write to gzip file " << r2_path << std::endl;
          return;
        }
    }
  }
};


class InputDistributor
{
public:
  InputDistributor(std::string path1, std::string path2, uint32_t total_reads)
  : gzr1_(path1), gzr2_(path2), total_reads_(total_reads) {}

  void inputThread1()
  {
    while (next_seq1_ < total_reads_)
    {
      uint32_t this_seq = next_seq1_++;

      std::vector<int> wanters;
      for (int i=0; i<the_consumers_.size(); i++)
        if (the_consumers_[i]->wantThisRead1(this_seq))
          wanters.push_back(i);

      if (wanters.empty())
        gzr1_.skip4lines();
      else
      {
        FastqHalfRead the_read(gzr1_);
        for (int i : wanters)
          the_consumers_[i]->q1.enqueue(the_read);
      }

      static int line_count = 0;
      line_count += 4;
      if (line_count % 4000000 > 3999995)
        std::cerr << "Processed " << (line_count / 4) << " reads..." << std::endl;
    }
    for (SubsampleWriter* c : the_consumers_)
      c->q1.shutdown();
  }
  void inputThread2()
  {
    while (next_seq2_ < total_reads_)
    {
      uint32_t this_seq = next_seq2_++;

      std::vector<int> wanters;
      for (int i=0; i<the_consumers_.size(); i++)
        if (the_consumers_[i]->wantThisRead2(this_seq))
          wanters.push_back(i);

      if (wanters.empty())
        gzr2_.skip4lines();
      else
      {
        FastqHalfRead the_read(gzr2_);
        for (int i : wanters)
          the_consumers_[i]->q2.enqueue(the_read);
      }
    }
    for (SubsampleWriter* c : the_consumers_)
      c->q2.shutdown();
  }

  void register_consumer(SubsampleWriter* w)
  {
    the_consumers_.push_back(w);
  }

  std::vector<SubsampleWriter*> the_consumers_;
  GzipReader gzr1_;
  GzipReader gzr2_;
  uint32_t total_reads_;
  uint32_t next_seq1_ = 0;
  uint32_t next_seq2_ = 0;
};

void randomSampleIndices(std::vector<uint32_t>& samples, uint32_t k, uint32_t n)
{
  samples.resize(n);
  for (uint32_t i = 0; i < n; i++)
    samples[i] = i;
  std::random_device rd;
  std::mt19937 gen(rd());
  for (uint32_t i = n - 1; i > 0; i--)
  {
    std::uniform_int_distribution<> dis(0, i);
    int j = dis(gen);
    std::swap(samples[i], samples[j]);
  }
  samples.resize(k);
  samples.shrink_to_fit();
  std::sort(samples.begin(), samples.end());
  std::cerr << "finished generating "<<k<<" random samples out of "<<n
            <<", first read index is "<<samples.front()
            <<" middle read index is "<<samples[samples.size()/2]
            <<" last read index is "<<samples.back()<<std::endl;
}

void crash(std::string msg) { std::cerr<<msg<<std::endl; exit(1); }

void allDigitsOrDie(std::string s)
{
  for (char c : s)
    if (c < '0' || c > '9')
      crash(std::string("Not a number: ") + s);
}

// parse CSV to determine outputs, and pick the random indices
std::vector<std::unique_ptr<SubsampleWriter>> configureWriters(std::string config_path,
                                                               uint32_t total_reads)
{
  std::ifstream config_file(config_path);
  if (!config_file)
    crash(std::string("Failed to open config file: ")+config_path);
  std::string line;
  std::vector<std::string> base_outpaths;
  std::vector<uint32_t> num_subsamples;
  while (std::getline(config_file, line))
  {
    std::stringstream ss(line);
    std::string base_path;
    std::string num_subsamples_str;
    if (std::getline(ss, base_path, ',') && std::getline(ss, num_subsamples_str))
    {
      while (!base_path.empty() && base_path.back() == ' ')
        base_path.pop_back();
      base_outpaths.push_back(base_path);
      while (!num_subsamples_str.empty() && num_subsamples_str.front() == ' ')
        num_subsamples_str = num_subsamples_str.substr(1);
      allDigitsOrDie(num_subsamples_str);
      num_subsamples.push_back(std::stoi(num_subsamples_str));
    }
    else
      crash(std::string("failed to parse csv line:\n")+line);
  }
  std::vector<std::unique_ptr<SubsampleWriter>> outputters(base_outpaths.size());
  for (int i=0; i<base_outpaths.size(); i++)
  {
    outputters[i] = std::make_unique<SubsampleWriter>(base_outpaths[i]+"r1.fastq.gz",
                                                      base_outpaths[i]+"r2.fastq.gz");
  }
  for (int i=0; i<outputters.size(); i++)
    randomSampleIndices(outputters[i]->selected_indices, num_subsamples[i], total_reads);

  return outputters;
}

int main(int argc, char* argv[])
{
  if (argc != 5)
  {
    std::cerr << "Usage: " << argv[0] << " <r1_file.fastq.gz> <r2_file.fastq.gz> "
              << "<total_reads_in_input_fastq> <config_csv>" << std::endl;
    return 1;
  }

  try
  {
    std::string r1_path = argv[1];
    std::string r2_path = argv[2];
    allDigitsOrDie(argv[3]);
    uint32_t total_reads = std::stoi(argv[3]);
    std::string config_path = argv[4];

    InputDistributor inputter(r1_path, r2_path, total_reads);
    auto outputters = configureWriters(config_path, total_reads);

    std::vector<std::thread> threads;
    for (auto& cur : outputters)
      inputter.register_consumer(cur.get());
    std::thread input_thread1([&inputter](){inputter.inputThread1();});
    std::thread input_thread2([&inputter](){inputter.inputThread2();});
    for (int i = 0; i < outputters.size(); i++)
    {
      SubsampleWriter* cur = outputters[i].get();
      threads.emplace_back([cur]() { cur->outputThread1(); });
      threads.emplace_back([cur]() { cur->outputThread2(); });
    }
    for (auto& t : threads)
      t.join();
    input_thread1.join();
    input_thread2.join();
  }
  catch (const std::exception& e) { crash(std::string("Error: ")+e.what()); }
  return 0;
}
