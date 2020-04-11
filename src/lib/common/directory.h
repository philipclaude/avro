#ifndef avro_LIB_COMMON_DIRECTORY_H_
#define avro_LIB_COMMON_DIRECTORY_H_

#include "common/error.h"
#include "common/json.h"

#include <json/json.hpp>

#include <dirent.h>
#include <map>
#include <string>
#include <vector>

#ifdef WINDOWS
  #include <direct.h>
  #define get_current_dir _getcwd
  #define filesep '\\'
#else
  #include <unistd.h>
  #define get_current_dir getcwd
  #define filesep '/'
#endif

inline std::string
get_file_ext( const std::string& filename )
{
  std::string::size_type idx;
  idx = filename.rfind('.'); // find the '.' in reverse order
  if (idx!=std::string::npos)
    return filename.substr(idx+1);
  return "";
}

namespace avro
{

class Directory
{

public:
  Directory( const std::string& initial=std::string() )
  {
    initialize(initial);
  }

  void initialize(const std::string& initial)
  {
    current_ = initial;
    if (current_.empty())
    {
      char current_path[FILENAME_MAX];
      if (!get_current_dir(current_path,sizeof(current_path)))
      {
        printf("error getting current directory");
        avro_assert_not_reached;
      }
      current_path[sizeof(current_path)-1] = '\0';
      current_.assign(current_path);
    }
  }

  void ls( std::vector<json>& directory , std::string path=std::string() )
  {
    if (path.empty()) path = current_;
    DIR *dir;
    struct dirent *ent;

    if ((dir=opendir(path.c_str()))!=NULL)
    {
      while ((ent=readdir(dir))!=NULL)
      {
        std::string s = std::string(ent->d_name);
        std::string type;
        if (ent->d_type==DT_DIR) type = "dir";
        else if (ent->d_type==DT_REG) type = "file";
        else continue;
        if (s==".") continue;
        json entry;
        entry["entry"] = s;
        entry["type"]  = type;
        directory.push_back(entry);
      }
      closedir(dir);
    }
    else
    {
      perror("");
    }
  }

  void cd( const std::string& directory )
  {
    if (directory=="..")
    {
      size_t slash = current_.rfind(filesep);
      current_ = current_.substr(0,slash);
    }
    else
      current_ += "/" + directory + "\0";
  }

  void clear()
  {
    items_.clear();
    current_ = "";
  }

  std::string pwd() const { return current_; }

private:
  enum DirectoryItemType { Folder , File , Path };

  std::map<std::string,DirectoryItemType> items_;
  std::string current_;
};

} // avro

#endif
