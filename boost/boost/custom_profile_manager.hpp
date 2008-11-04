
// This is a demo of a custom profiler 

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include <profiler/profiler.hpp>

using namespace boost::multi_index;
using namespace boost::prof;

struct meta_profile_data 
{
  meta_profile_data(string n, string f, int l) : name(n), file(f), line(l) { };
  meta_profile_data(const meta_profile_data& x) : name(x.name), file(x.file), line(x.line) { };
  bool operator<(const meta_profile_data& x) const { 
    return name < x.name; 
  }
  string name;
  string file;
  int line;
  
  friend ostream& operator<<(ostream& o, const meta_profile_data& x)  { 
    return o << x.name << ":" << x.file << ":" << x.line; 
  }
};


typedef basic_profile_manager
  <
    meta_profile_data,
    double, 
    boost::high_resolution_timer,
    default_logging_policy,
    default_stats_policy<meta_profile_data, double>
  > 
  custom_profile_manager;

typedef custom_profile_manager::stats_policy::stats_map profile_stats;
typedef basic_profiler<custom_profile_manager> custom_profiler;

profile_stats& GetProfileStats() {
  return custom_profile_manager::stats_policy::stats; 
}

struct meta_stats {
  meta_stats(const pair<meta_profile_data, counted_sum<double> >& x)
    : name(x.first.name)
    , file(x.first.file)
    , line(x.first.line)
    , count(x.second.first)
    , dur(x.second.second) 
  { }
  string name;
  string file;
  int line;
  int count;
  double dur;
};

ostream& operator<<(ostream& o, const meta_stats& x) { 
  return o << x.name << " " /*<< x.file << " "*/ << x.line << " " << x.count << " " << x.dur << endl; 
}

struct index_name { };
struct index_dur { };

typedef multi_index_container<
  meta_stats,
  indexed_by<
    ordered_non_unique<tag<index_name>,  member< meta_stats, string, &meta_stats::name> >,
    ordered_non_unique<tag<index_dur>, member< meta_stats, double, &meta_stats::dur> >
  >
> meta_stats_container;

template<typename index_t>
void OutputOrderedMetaStats(const meta_stats_container& x) {
  const typename boost::multi_index::index<meta_stats_container,index_t>::type& i = 
    boost::multi_index::get<index_t>(x);
  std::copy(i.begin(),i.end(),std::ostream_iterator<meta_stats_container::value_type>(std::cout));
}

void GenerateCustomReport() 
{
  meta_stats_container cont; 
  profile_stats::const_iterator i = GetProfileStats().begin();
  while (i != GetProfileStats().end()) {
    cont.insert(meta_stats(*i++));
  }

  cout << endl << "Stats ordered by name " << endl;
  OutputOrderedMetaStats<index_name>(cont);

  cout << endl << "Stats ordered by duration " << endl;
  OutputOrderedMetaStats<index_dur>(cont);
}
