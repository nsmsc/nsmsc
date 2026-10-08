#include "nsmsc/rr_set_oracle.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <fstream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>

namespace nsmsc {
namespace {
constexpr std::array<char, 8> kMagic{'N','S','R','R','S','E','T','1'};
constexpr std::array<char, 8> kComponentMagic{'N','S','R','R','C','M','P','1'};
template<class T> void writeLe(std::ostream& out, T value) {
    static_assert(std::is_unsigned_v<T>);
    for (std::size_t i=0;i<sizeof(T);++i) out.put(static_cast<char>((value>>(8*i))&0xff));
}
template<class T> T readLe(std::istream& in) {
    static_assert(std::is_unsigned_v<T>); T value=0;
    for (std::size_t i=0;i<sizeof(T);++i) { int c=in.get(); if(c==EOF) throw std::runtime_error("Truncated RR file"); value|=static_cast<T>(static_cast<unsigned char>(c))<<(8*i); }
    return value;
}
}

InfluenceGraph readHiggsRetweetGraph(const std::filesystem::path& path) {
    std::ifstream input(path); if (!input) throw std::runtime_error("Cannot open " + path.string());
    struct Edge { std::uint64_t source, actor; };
    std::vector<Edge> edges; std::vector<std::uint64_t> ids;
    std::uint64_t actor, source, weight;
    while (input >> actor >> source >> weight) {
        if (actor == source) continue;
        edges.push_back({source, actor}); ids.push_back(source); ids.push_back(actor);
    }
    if (!input.eof()) throw std::runtime_error("Malformed Higgs retweet graph");
    std::sort(ids.begin(),ids.end()); ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
    if (ids.size() > std::numeric_limits<ElementId>::max()) throw std::runtime_error("Too many graph nodes");
    std::unordered_map<std::uint64_t,ElementId> index; index.reserve(ids.size());
    for (std::size_t i=0;i<ids.size();++i) index.emplace(ids[i],static_cast<ElementId>(i));
    InfluenceGraph graph{ids,std::vector<std::vector<ElementId>>(ids.size())};
    for (const auto& edge:edges) graph.incoming[index[edge.actor]].push_back(index[edge.source]);
    for (auto& list:graph.incoming) { std::sort(list.begin(),list.end()); list.erase(std::unique(list.begin(),list.end()),list.end()); }
    return graph;
}

WeightedCascadeSampler::WeightedCascadeSampler(
    const InfluenceGraph& graph,
    std::uint64_t seed,
    bool subsimJumpSampling)
    : graph_(graph), random_(seed), marks_(graph.externalIds.size()),
      subsimJumpSampling_(subsimJumpSampling) {}

void WeightedCascadeSampler::append(RRSetCollection& collection, std::size_t count) {
    if (collection.sets.empty()) { collection.nodeCount=graph_.externalIds.size(); }
    if (collection.nodeCount!=graph_.externalIds.size()) throw std::invalid_argument("RR graph size mismatch");
    std::uniform_int_distribution<std::size_t> root(0,collection.nodeCount-1);
    std::uniform_real_distribution<double> unit(0.0,1.0);
    for (std::size_t sample=0;sample<count;++sample) {
        if (++generation_==0) { std::fill(marks_.begin(),marks_.end(),0); ++generation_; }
        std::vector<ElementId> rr; std::queue<ElementId> queue;
        const auto start=static_cast<ElementId>(root(random_)); marks_[start]=generation_; queue.push(start);
        while(!queue.empty()) {
            const ElementId v=queue.front(); queue.pop(); rr.push_back(v);
            const auto& incoming=graph_.incoming[v];
            const double probability=incoming.empty()?0.0:1.0/static_cast<double>(incoming.size());
            if (subsimJumpSampling_ && probability < 1.0 && !incoming.empty()) {
                // SUBSIM's WC specialization samples successful edge positions
                // geometrically instead of tossing one coin for every edge.
                std::geometric_distribution<std::size_t> skip(probability);
                std::size_t position = skip(random_);
                while (position < incoming.size()) {
                    const ElementId u = incoming[position];
                    if (marks_[u] != generation_) {
                        marks_[u] = generation_;
                        queue.push(u);
                    }
                    position += 1 + skip(random_);
                }
            } else {
                for (ElementId u:incoming) if(marks_[u]!=generation_ && unit(random_)<probability) { marks_[u]=generation_; queue.push(u); }
            }
        }
        std::sort(rr.begin(),rr.end()); collection.sets.push_back(std::move(rr));
    }
}

void writeRRSetsBinary(const std::filesystem::path& path,const RRSetCollection& c) {
    std::ofstream out(path,std::ios::binary|std::ios::trunc); if(!out) throw std::runtime_error("Cannot open " + path.string());
    out.write(kMagic.data(),kMagic.size()); writeLe(out,std::uint64_t(c.nodeCount)); writeLe(out,c.seed); writeLe(out,std::uint64_t(c.sets.size()));
    for(const auto& set:c.sets) { writeLe(out,std::uint64_t(set.size())); for(ElementId e:set) writeLe(out,e); }
    if(!out) throw std::runtime_error("Failed to write " + path.string());
}
RRSetCollection readRRSetsBinary(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary); if(!in) throw std::runtime_error("Cannot open " + path.string());
    std::array<char,8> magic{}; in.read(magic.data(),magic.size()); if(magic!=kMagic) throw std::runtime_error("Invalid RR file");
    RRSetCollection c; c.nodeCount=readLe<std::uint64_t>(in); c.seed=readLe<std::uint64_t>(in); const auto count=readLe<std::uint64_t>(in); c.sets.reserve(count);
    for(std::uint64_t i=0;i<count;++i) { const auto size=readLe<std::uint64_t>(in); std::vector<ElementId> set(size); for(auto& e:set) { e=readLe<ElementId>(in); if(e>=c.nodeCount) throw std::runtime_error("RR element out of range"); } c.sets.push_back(std::move(set)); }
    if(in.peek()!=EOF) throw std::runtime_error("RR file has trailing bytes");
    return c;
}

void writeRRComplementarityBinary(const std::filesystem::path& path,const RRComplementarityConfig& c) {
    if(c.nodeCount==0||c.tau==0||c.rho<0)throw std::invalid_argument("Invalid RR complementarity config");
    std::ofstream out(path,std::ios::binary|std::ios::trunc);if(!out)throw std::runtime_error("Cannot open "+path.string());
    out.write(kComponentMagic.data(),kComponentMagic.size());writeLe(out,std::uint64_t(c.nodeCount));writeLe(out,std::uint64_t(c.tau));writeLe(out,std::bit_cast<std::uint64_t>(c.rho));writeLe(out,std::uint64_t(c.components.size()));
    for(const auto& component:c.components){if(component.witness.empty()||component.witness.size()>c.tau||component.value<0)throw std::invalid_argument("Invalid RR component");writeLe(out,std::bit_cast<std::uint64_t>(component.value));writeLe(out,std::uint64_t(component.witness.size()));for(ElementId e:component.witness){if(e>=c.nodeCount)throw std::invalid_argument("RR witness out of range");writeLe(out,e);}}
    if(!out)throw std::runtime_error("Failed to write "+path.string());
}
RRComplementarityConfig readRRComplementarityBinary(const std::filesystem::path& path,std::size_t expectedNodeCount){
    std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("Cannot open "+path.string());std::array<char,8> magic{};in.read(magic.data(),magic.size());if(magic!=kComponentMagic)throw std::runtime_error("Invalid RR component file");
    RRComplementarityConfig c;c.nodeCount=readLe<std::uint64_t>(in);c.tau=readLe<std::uint64_t>(in);c.rho=std::bit_cast<double>(readLe<std::uint64_t>(in));const auto h=readLe<std::uint64_t>(in);if(c.nodeCount!=expectedNodeCount||c.tau==0||c.rho<0)throw std::runtime_error("RR component metadata mismatch");c.components.reserve(h);
    for(std::uint64_t i=0;i<h;++i){ComplementaryComponent component;component.value=std::bit_cast<double>(readLe<std::uint64_t>(in));const auto size=readLe<std::uint64_t>(in);if(size==0||size>c.tau)throw std::runtime_error("Invalid RR witness size");component.witness.resize(size);for(auto& e:component.witness){e=readLe<ElementId>(in);if(e>=c.nodeCount)throw std::runtime_error("RR witness out of range");}c.components.push_back(std::move(component));}if(in.peek()!=EOF)throw std::runtime_error("RR component file has trailing bytes");return c;
}

RRSetOracle::RRSetOracle(RRSetCollection c,std::vector<ComplementaryComponent> components)
    : collection_(std::move(c)),components_(std::move(components)),setsByElement_(collection_.nodeCount),componentsByElement_(collection_.nodeCount),coverCount_(collection_.sets.size()),selected_(collection_.nodeCount) {
    if(collection_.sets.empty()) throw std::invalid_argument("RR collection must not be empty");
    for(std::size_t i=0;i<collection_.sets.size();++i) for(ElementId e:collection_.sets[i]) { validateElement(e); setsByElement_[e].push_back(static_cast<std::uint32_t>(i)); }
    for(std::size_t i=0;i<components_.size();++i) {
        const auto& cpt=components_[i];
        for(std::size_t j=0;j<cpt.witness.size();++j) {
            const ElementId e=cpt.witness[j];
            validateElement(e);
            if(std::find(cpt.witness.begin(),cpt.witness.begin()+j,e)==cpt.witness.begin()+j)
                componentsByElement_[e].push_back(static_cast<ComponentId>(i));
        }
    }
}
std::size_t RRSetOracle::groundSize() const{return selected_.size();} std::size_t RRSetOracle::componentCount() const{return components_.size();}
void RRSetOracle::validateElement(ElementId e) const{if(e>=selected_.size()) throw std::out_of_range("Element out of range");}
void RRSetOracle::reset(){std::fill(selected_.begin(),selected_.end(),0);std::fill(coverCount_.begin(),coverCount_.end(),0);selectedSize_=coveredSets_=0;}
bool RRSetOracle::selected(ElementId e) const{validateElement(e);return selected_[e];} std::size_t RRSetOracle::selectedSize() const{return selectedSize_;}
double RRSetOracle::fValue() const{return static_cast<double>(groundSize())*coveredSets_/collection_.sets.size();}
double RRSetOracle::fMarginal(ElementId e) const{validateElement(e);if(selected_[e])return 0;std::size_t gain=0;for(auto r:setsByElement_[e])gain+=coverCount_[r]==0;return static_cast<double>(groundSize())*gain/collection_.sets.size();}
double RRSetOracle::totalMarginal(ElementId e) const{if(selected(e))return 0;double gain=fMarginal(e);for(ComponentId i:componentsByElement_[e]){const auto& c=components_[i];bool completes=true;for(ElementId w:c.witness)if(w!=e&&!selected_[w]){completes=false;break;}if(completes)gain+=c.value;}return gain;}
void RRSetOracle::add(ElementId e){validateElement(e);if(selected_[e])return;selected_[e]=1;++selectedSize_;for(auto r:setsByElement_[e])if(coverCount_[r]++==0)++coveredSets_;}
void RRSetOracle::remove(ElementId e){validateElement(e);if(!selected_[e])return;selected_[e]=0;--selectedSize_;for(auto r:setsByElement_[e])if(--coverCount_[r]==0)--coveredSets_;}
double RRSetOracle::gValue(ComponentId i) const{if(i>=components_.size())throw std::out_of_range("Component out of range");for(ElementId e:components_[i].witness)if(!selected_[e])return 0;return components_[i].value;}
double RRSetOracle::gMaximum(ComponentId i) const{if(i>=components_.size())throw std::out_of_range("Component out of range");return components_[i].value;}
std::span<const ElementId> RRSetOracle::witness(ComponentId i) const{if(i>=components_.size())throw std::out_of_range("Component out of range");return components_[i].witness;}
double RRSetOracle::totalValue() const{double result=fValue();for(std::size_t i=0;i<components_.size();++i)result+=gValue(static_cast<ComponentId>(i));return result;}
}
