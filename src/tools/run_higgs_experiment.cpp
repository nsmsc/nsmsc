#include "nsmsc/baselines.hpp"
#include "nsmsc/complementary_completion.hpp"
#include "nsmsc/experiment_protocol.hpp"
#include "nsmsc/greedy.hpp"
#include "nsmsc/rr_set_oracle.hpp"

#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace {
std::uint64_t number(std::string_view s,const char* name){std::uint64_t v=0;auto r=std::from_chars(s.data(),s.data()+s.size(),v);if(r.ec!=std::errc{}||r.ptr!=s.data()+s.size())throw std::invalid_argument(std::string("Invalid ")+name);return v;}
double realNumber(std::string_view text, const char* name) {
    std::size_t consumed = 0;
    const std::string value(text);
    const double result = std::stod(value, &consumed);
    if (consumed != value.size()) throw std::invalid_argument(std::string("Invalid ") + name);
    return result;
}
bool methodEnabled(
    std::string_view filter, std::string_view method) {
    if (filter.empty()) return true;
    std::size_t begin = 0;
    while (begin <= filter.size()) {
        const std::size_t end = filter.find(',', begin);
        const auto token = filter.substr(
            begin, end == std::string_view::npos ? end : end - begin);
        if (token == method) return true;
        if (end == std::string_view::npos) break;
        begin = end + 1;
    }
    return false;
}
struct Result{std::string method,group;std::size_t budget,augmentation;nsmsc::Solution solution;};
struct Evaluation{double f=0,g=0,total=0;};
Evaluation evaluate(const nsmsc::RRSetCollection& pool,const std::vector<nsmsc::ComplementaryComponent>& components,const std::vector<nsmsc::ElementId>& elements){std::vector<std::uint8_t> selected(pool.nodeCount);for(auto e:elements)selected[e]=1;std::size_t covered=0;for(const auto& set:pool.sets)covered+=std::any_of(set.begin(),set.end(),[&](auto e){return selected[e]!=0;});const double f=static_cast<double>(pool.nodeCount)*covered/pool.sets.size();double g=0;for(const auto& component:components)if(std::all_of(component.witness.begin(),component.witness.end(),[&](auto e){return selected[e]!=0;}))g+=component.value;return {f,g,f+g};}

std::vector<nsmsc::ElementId> readOrder(
    const std::filesystem::path& path, std::size_t nodeCount) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open seed order " + path.string());
    std::vector<nsmsc::ElementId> order;
    std::unordered_set<nsmsc::ElementId> seen;
    std::uint64_t value = 0;
    while (input >> value) {
        if (value >= nodeCount ||
            !seen.insert(static_cast<nsmsc::ElementId>(value)).second) {
            throw std::runtime_error("Invalid seed order " + path.string());
        }
        order.push_back(static_cast<nsmsc::ElementId>(value));
    }
    if (!input.eof()) throw std::runtime_error("Malformed seed order " + path.string());
    return order;
}

double readAuthorRuntime(const std::filesystem::path& seedPath) {
    const std::string seedName = seedPath.filename().string();
    if (!seedName.starts_with("seed_")) return 0.0;
    const auto resultPath = seedPath.parent_path().parent_path() /
        seedName.substr(5);
    std::ifstream input(resultPath);
    std::string label;
    double runtime = 0.0;
    while (input >> label) {
        if (label == "Time") {
            std::string seconds;
            input >> seconds >> runtime;
            return runtime;
        }
        std::string rest;
        std::getline(input, rest);
    }
    return 0.0;
}

nsmsc::Solution solutionFromOrder(
    nsmsc::ObjectiveOracle& oracle,
    const std::vector<nsmsc::ElementId>& order,
    std::size_t budget,
    double runtime) {
    if (order.size() < budget) throw std::invalid_argument("Seed order is too short");
    oracle.reset();
    nsmsc::Solution solution;
    solution.elements.assign(order.begin(), order.begin() + budget);
    for (auto element : solution.elements) oracle.add(element);
    solution.fValue = oracle.fValue();
    solution.complementaryValue = oracle.totalValue() - solution.fValue;
    solution.totalValue = oracle.totalValue();
    solution.runtimeSeconds = runtime;
    return solution;
}

struct ExternalOrders {
    std::filesystem::path histK, histAugmented, opimcK, opimcAugmented;
};
}

int main(int argc,char** argv){
    if(argc<7){std::cerr<<"Usage: run_higgs_experiment <optimization.rrbin> <evaluation.rrbin> <components.bin> <k> <output.jsonl> <algorithm-seed> [--dataset=name] [--augmentation=b] [--hist-k=file] [--hist-aug=file] [--opimc-k=file] [--opimc-aug=file] [--ls-max-evals=n] [--ls-max-seconds=s]\n";return 2;}
    try{
        auto optimization=nsmsc::readRRSetsBinary(argv[1]);const auto evaluation=nsmsc::readRRSetsBinary(argv[2]);if(optimization.nodeCount!=evaluation.nodeCount)throw std::invalid_argument("RR pool node counts differ");
        const auto config=nsmsc::readRRComplementarityBinary(argv[3],optimization.nodeCount);const auto k=static_cast<std::size_t>(number(argv[4],"k"));const auto seed=number(argv[6],"algorithm seed");const auto theoretical=nsmsc::theoreticalAugmentation(config.components.size(),config.tau);std::size_t augmentation=theoretical;ExternalOrders external;std::uint64_t lsMaxEvaluations=10'000'000;double lsMaxSeconds=60.0;std::string datasetName="Higgs";
        for(int i=7;i<argc;++i){const std::string_view option(argv[i]);const auto value=[&](std::string_view prefix){return std::string(option.substr(prefix.size()));};if(option.starts_with("--dataset="))datasetName=value("--dataset=");else if(option.starts_with("--augmentation="))augmentation=number(option.substr(15),"augmentation");else if(option.starts_with("--hist-k="))external.histK=value("--hist-k=");else if(option.starts_with("--hist-aug="))external.histAugmented=value("--hist-aug=");else if(option.starts_with("--opimc-k="))external.opimcK=value("--opimc-k=");else if(option.starts_with("--opimc-aug="))external.opimcAugmented=value("--opimc-aug=");else if(option.starts_with("--ls-max-evals="))lsMaxEvaluations=number(option.substr(15),"LS max evaluations");else if(option.starts_with("--ls-max-seconds="))lsMaxSeconds=realNumber(option.substr(17),"LS max seconds");else throw std::invalid_argument("Unknown option: "+std::string(option));}if(datasetName.empty())throw std::invalid_argument("Dataset name must not be empty");if(k+augmentation>optimization.nodeCount)throw std::invalid_argument("Budget exceeds node count");const auto augmented=k+augmentation;
        const auto optimizationSetCount=optimization.sets.size();const char* filterValue=std::getenv("METHOD_FILTER");const std::string filter=filterValue?filterValue:"";const auto enabled=[&](std::string_view method){return methodEnabled(filter,method);};nsmsc::RRSetOracle oracle(std::move(optimization),config.components);std::vector<Result> results;
        const auto pair=[&](const std::string& method,auto run){if(!enabled(method))return;results.push_back({method,"augmented",augmented,augmentation,run(augmented)});results.push_back({method,"original",k,0,run(k)});};
        pair("FTG-2022",[&](std::size_t b){return nsmsc::fastThresholdGreedyOnTotal(oracle,b,0.1,1.0);});
        pair("F-Greedy-LS",[&](std::size_t b){return nsmsc::fGreedyLocalSearch(oracle,b,0.1);});
        pair("RRG",[&](std::size_t b){return nsmsc::residualRandomGreedy(oracle,b,seed);});
        pair("threshold_greedy_2024",[&](std::size_t b){return nsmsc::thresholdGreedyOnTotal(oracle,b,0.1);});
        if(enabled("proposed")){results.push_back({"proposed","augmented",augmented,augmentation,nsmsc::greedyCompletion(oracle,k,augmentation,config.tau)});results.push_back({"proposed","original",k,0,nsmsc::greedyCompletion(oracle,k,0,config.tau)});}
        if(enabled("proposed-original")){results.push_back({"proposed-original","augmented",augmented,augmentation,nsmsc::greedyCompletionOriginal(oracle,k,augmentation,config.tau)});results.push_back({"proposed-original","original",k,0,nsmsc::greedyCompletionOriginal(oracle,k,0,config.tau)});}
        if (enabled("proposed-no-enumeration")) {
            results.push_back({"proposed-no-enumeration", "augmented", augmented,
                augmentation, nsmsc::greedyCompletionFixedQ(
                    oracle, k, augmentation, config.tau)});
            results.push_back({"proposed-no-enumeration", "original", k, 0,
                nsmsc::greedyCompletionFixedQ(oracle, k, 0, config.tau)});
        }
        if (enabled("proposed-no-ranking")) {
            results.push_back({"proposed-no-ranking", "augmented", augmented,
                augmentation, nsmsc::greedyCompletionRandomRanking(
                    oracle, k, augmentation, config.tau, seed)});
            results.push_back({"proposed-no-ranking", "original", k, 0,
                nsmsc::greedyCompletionRandomRanking(
                    oracle, k, 0, config.tau, seed)});
        }
        if(enabled("proposed-threshold")){results.push_back({"proposed-threshold","augmented",augmented,augmentation,nsmsc::thresholdTargetedCompletion(oracle,k,augmentation,config.tau)});results.push_back({"proposed-threshold","original",k,0,nsmsc::thresholdTargetedCompletion(oracle,k,0,config.tau)});}
        if(enabled("proposed-cost-aware")){results.push_back({"proposed-cost-aware","augmented",augmented,augmentation,nsmsc::costAwareCompletion(oracle,k,augmentation,config.tau)});results.push_back({"proposed-cost-aware","original",k,0,nsmsc::costAwareCompletion(oracle,k,0,config.tau)});}
        if(enabled("proposed-adaptive")){auto solution=nsmsc::adaptiveCompletion(oracle,k,augmentation,config.tau);const auto used=solution.elements.size()-k;results.push_back({"proposed-adaptive","augmented",k+used,used,std::move(solution)});results.push_back({"proposed-adaptive","original",k,0,nsmsc::adaptiveCompletion(oracle,k,0,config.tau)});}
        if (enabled("HIST") || enabled("proposed-HIST") ||
            enabled("proposed-adaptive-HIST")) {
            if (external.histK.empty() || external.histAugmented.empty())
                throw std::invalid_argument("HIST seed files are required");
            const auto histK = readOrder(external.histK, oracle.groundSize());
            const auto histAug = readOrder(external.histAugmented, oracle.groundSize());
            const double timeK = readAuthorRuntime(external.histK);
            const double timeAug = readAuthorRuntime(external.histAugmented);
            if (enabled("HIST")) {
                results.push_back({"HIST", "augmented", augmented, augmentation,
                    solutionFromOrder(oracle, histAug, augmented, timeAug)});
                results.push_back({"HIST", "original", k, 0,
                    solutionFromOrder(oracle, histK, k, timeK)});
            }
            if (enabled("proposed-HIST")) {
                auto augResult = nsmsc::greedyCompletionFromFOrder(
                    oracle, k, augmentation, config.tau, histAug);
                augResult.runtimeSeconds += timeAug;
                results.push_back({"proposed-HIST", "augmented", augmented,
                    augmentation, std::move(augResult)});
                auto baseResult = nsmsc::greedyCompletionFromFOrder(
                    oracle, k, 0, config.tau, histK);
                baseResult.runtimeSeconds += timeK;
                results.push_back({"proposed-HIST", "original", k, 0,
                    std::move(baseResult)});
            }
            if (enabled("proposed-adaptive-HIST")) {
                auto adapted = nsmsc::adaptiveCompletionFromFOrder(
                    oracle, k, augmentation, config.tau, histK, histAug);
                adapted.runtimeSeconds += timeK + timeAug;
                const auto used = adapted.elements.size() - k;
                results.push_back({"proposed-adaptive-HIST", "augmented",
                    k + used, used, std::move(adapted)});
                auto original = nsmsc::adaptiveCompletionFromFOrder(
                    oracle, k, 0, config.tau, histK, histK);
                original.runtimeSeconds += timeK;
                results.push_back({"proposed-adaptive-HIST", "original",
                    k, 0, std::move(original)});
            }
        }
        if(enabled("HIST-LS")||enabled("proposed-threshold-HIST")||enabled("proposed-cost-aware-HIST")){if(external.histK.empty()||external.histAugmented.empty())throw std::invalid_argument("HIST seed files are required");const auto histK=readOrder(external.histK,oracle.groundSize());const auto histAug=readOrder(external.histAugmented,oracle.groundSize());const double timeK=readAuthorRuntime(external.histK),timeAug=readAuthorRuntime(external.histAugmented);if(enabled("HIST-LS")){auto augLs=nsmsc::localSearchFromOrder(oracle,histAug,augmented,0.1,lsMaxEvaluations,lsMaxSeconds);augLs.runtimeSeconds+=timeAug;results.push_back({"HIST-LS","augmented",augmented,augmentation,std::move(augLs)});auto baseLs=nsmsc::localSearchFromOrder(oracle,histK,k,0.1,lsMaxEvaluations,lsMaxSeconds);baseLs.runtimeSeconds+=timeK;results.push_back({"HIST-LS","original",k,0,std::move(baseLs)});}if(enabled("proposed-threshold-HIST")){auto augResult=nsmsc::thresholdTargetedCompletionFromFOrder(oracle,k,augmentation,config.tau,histK);augResult.runtimeSeconds+=timeK;results.push_back({"proposed-threshold-HIST","augmented",augmented,augmentation,std::move(augResult)});auto baseResult=nsmsc::thresholdTargetedCompletionFromFOrder(oracle,k,0,config.tau,histK);baseResult.runtimeSeconds+=timeK;results.push_back({"proposed-threshold-HIST","original",k,0,std::move(baseResult)});}if(enabled("proposed-cost-aware-HIST")){auto augResult=nsmsc::costAwareCompletionFromFOrder(oracle,k,augmentation,config.tau,histK);augResult.runtimeSeconds+=timeK;results.push_back({"proposed-cost-aware-HIST","augmented",augmented,augmentation,std::move(augResult)});auto baseResult=nsmsc::costAwareCompletionFromFOrder(oracle,k,0,config.tau,histK);baseResult.runtimeSeconds+=timeK;results.push_back({"proposed-cost-aware-HIST","original",k,0,std::move(baseResult)});}}
        if(enabled("OPIM-C")){if(external.opimcK.empty()||external.opimcAugmented.empty())throw std::invalid_argument("OPIM-C seed files are required");const auto opimK=readOrder(external.opimcK,oracle.groundSize());const auto opimAug=readOrder(external.opimcAugmented,oracle.groundSize());results.push_back({"OPIM-C","augmented",augmented,augmentation,solutionFromOrder(oracle,opimAug,augmented,readAuthorRuntime(external.opimcAugmented))});results.push_back({"OPIM-C","original",k,0,solutionFromOrder(oracle,opimK,k,readAuthorRuntime(external.opimcK))});}
        pair("f_greedy",[&](std::size_t b){return nsmsc::greedyOnF(oracle,b);});pair("F_greedy",[&](std::size_t b){return nsmsc::greedyOnTotal(oracle,b);});pair("random",[&](std::size_t b){return nsmsc::randomSelection(oracle,b,seed);});
        if(results.empty())throw std::invalid_argument("Unknown METHOD_FILTER");const std::filesystem::path outputPath(argv[5]);if(outputPath.has_parent_path())std::filesystem::create_directories(outputPath.parent_path());std::ofstream output(outputPath);if(!output)throw std::runtime_error("Cannot create result file");
        for(const auto& result:results){const auto measured=evaluate(evaluation,config.components,result.solution.elements);output<<std::fixed<<std::setprecision(9)<<"{\"dataset\":\""<<datasetName<<"\",\"method\":\""<<result.method<<"\",\"comparison_group\":\""<<result.group<<"\",\"k\":"<<k<<",\"h\":"<<config.components.size()<<",\"tau\":"<<config.tau<<",\"rho\":"<<config.rho<<",\"augmentation\":"<<result.augmentation<<",\"allowed_budget\":"<<result.budget<<",\"selected_size\":"<<result.solution.elements.size()<<",\"optimization_f_value\":"<<result.solution.fValue<<",\"optimization_total_value\":"<<result.solution.totalValue<<",\"evaluation_f_value\":"<<measured.f<<",\"evaluation_complementary_value\":"<<measured.g<<",\"evaluation_total_value\":"<<measured.total<<",\"optimization_seconds\":"<<result.solution.runtimeSeconds<<",\"local_search_swaps\":"<<result.solution.localSearchSwaps<<",\"local_search_evaluations\":"<<result.solution.localSearchEvaluations<<",\"optimization_rr_sets\":"<<optimizationSetCount<<",\"evaluation_rr_sets\":"<<evaluation.sets.size()<<",\"algorithm_seed\":"<<seed<<",\"winning_branch\":\""<<result.solution.winningBranch<<"\"}\n";}
        if(!output)throw std::runtime_error("Failed to write results");std::cout<<"output="<<outputPath<<" runs="<<results.size()<<" h="<<config.components.size()<<" tau="<<config.tau<<" augmentation="<<augmentation<<'\n';
    }catch(const std::exception& e){std::cerr<<"run_higgs_experiment: "<<e.what()<<'\n';return 1;}return 0;
}
