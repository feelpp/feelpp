/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
   Date: 2023-07-24

   Copyright (C) 2023-present Feel++ Consortium
   

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#pragma once

#include <feel/feelcore/enumerate.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/table.hpp>
#include <feel/feelmor/crbdata.hpp>
#include <feel/feelmor/crbmodeldb.hpp>
#include <feel/feelmor/crbplugin_interface.hpp>
#include <fmt/ranges.h>

namespace Feel
{

/**
 * @brief Model wrapper for RB online
 *
 */
class MORModel
{
    public:
    std::string dbroot;
    fs::path dbroot_path;
    std::string name;
    std::string attribute = "last_modified"; // last_created, last_modified
    std::string load = "rb";
    std::string output;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE( MORModel, dbroot, name, attribute, load, output )

    template <typename... Args>
    auto
    run( Args&&... args ) const
    {
        return p_->run( std::forward<Args>( args )... );
    }
    auto
    parameterSpace() const
    {
        return p_->parameterSpace();
    }

    void
    initExporter()
    {
        if ( load == "all" && load == "fe" )
            p_->initExporter();
    }
    template <typename... Args>
    void
    exportField( Args&&... args ) const
    {
        if ( load == "all" || load == "fe" )
            p_->exportField( std::forward<Args>( args )... );
    }
    void
    saveExporter()
    {
        if ( load == "all" || load == "fe" )
            p_->saveExporter();
    }
    void
    loadPlugin();
    
private:
    std::shared_ptr<Feel::CRBPluginAPI> p_;
};

/**
 * @brief Observer for MOR online models
 *
 */
class MORObserver
{
  public:
    virtual ~MORObserver() = default;
    virtual void update( std::pair<int, ParameterSpaceX::Element> const&, std::vector<CRBResults> const& results ) = 0;
};
/**
 * @brief List of models for a given problem
 *
 * The models should at least share the same parameter space
 *
 */
class MORModels : public std::vector<MORModel>
{
    using super = std::vector<MORModel>;
    public:

    MORModels() = default;
    MORModels( nl::json && j ) : super( j["mormodels"].get<std::vector<MORModel>>() ) {}
    MORModels( nl::json const& j ) : super( j["mormodels"].get<std::vector<MORModel>>() ) {}

    void load()
    {
        for ( auto& m : *this )
        {
            m.loadPlugin();
        }
    }
    auto
    parameterSpace() const
    {
        return this->at( 0 ).parameterSpace();
    }

    /**
     * @brief a sampling of the parameter space
     *
     * @return auto
     */
    auto sampling() const
    {
        return this->at( 0 ).parameterSpace()->sampling();
    }
    /**
     * @brief initialize the Exporter for all models using the first model
     *
     */
    void
    initExporter()
    {
        if ( this->at( 0 ).load == "all" || this->at( 0 ).load == "fe" )
            this->at( 0 ).initExporter();
    }
    void
    saveExporter()
    {
        if ( this->at( 0 ).load == "all" || this->at( 0 ).load == "fe" )
            this->at( 0 ).saveExporter();
    }
    std::vector<std::vector<CRBResults>>
    run( std::shared_ptr<ParameterSpaceX::Sampling> const& sampling, nl::json const& data = {} ) const;

    void addObserver( std::shared_ptr<MORObserver> const& o )
    {
        observers_.push_back( o );
    }
private:    
    std::vector<std::shared_ptr<MORObserver>> observers_;
};

class MORExporter : public MORObserver
{
  public:
    MORExporter( MORModels& models )
        : models_( models )
    {
        models_.initExporter();
    }
    void update( std::pair<int, ParameterSpaceX::Element> const& mu, std::vector<CRBResults> const& results )
    {
        for ( auto const& [index, o] : enumerate( results[0].outputs() ) )
        {
            for ( auto const& p : models_ )
            {
                if ( p.load == "all" || p.load == "fe" )
                    p.exportField( fmt::format( "sol-{}-{}-{}", p.name, p.output, mu.first ), results[index] );
            }
        }
        models_.saveExporter();
    }
    MORModels& models_;
};
class MORTable : public MORObserver
{
  public:
    MORTable( MORModels const& models, bool print_to_screen = false, bool save_to_file = true )
        : models_( models ), print_to_screen_( print_to_screen ), save_to_file_( save_to_file ), results_path_( "" )
    {
        table_.format().setFloatingPointPrecision( ioption( _name = "output_results.precision" ) );
        std::vector<std::string> tableRowHeader; // = muspace->parameterNames();

        for ( auto const& p : models_ )
        {
            tableRowHeader.push_back( fmt::format( "{}.{} time", p.name, p.output ) );
            tableRowHeader.push_back( fmt::format( "{}.{} output", p.name, p.output ) );
            tableRowHeader.push_back( fmt::format( "{}.{} error", p.name, p.output ) );
        }
        LOG( INFO ) << fmt::format( "[MORModels::run] tableRowHeader {}", tableRowHeader );
        table_.add_row( tableRowHeader );
        table_.format().setFirstRowIsHeader( true );

        results_path_ = fs::path( models_.at( 0 ).dbroot_path ).parent_path() / fs::path( std::string( "results" ) );
        fs::create_directories( results_path_ );
        LOG( INFO ) << fmt::format( "[MORTable] results path : {}", results_path_.string() );
    }

    void update( std::pair<int, ParameterSpaceX::Element> const& mu, std::vector<CRBResults> const& results )
    {
        std::ostringstream ostrmu;
        for ( uint16_type d = 0; d < mu.second.size(); ++d )
            ostrmu << mu.second( d ) << " ";
        std::vector<double> tableRowValues( 3 * models_.size() );

        for ( auto const& [index, o] : enumerate( results[0].outputs() ) )
        {
            int curRowValIndex = 0;
            for ( auto const& [output_index, p] : enumerate( models_ ) )
            {
                int tableindex = output_index * 3;
                tableRowValues[tableindex + 0] = index;
                tableRowValues[tableindex + 1] = results[output_index].outputs()[index];
                tableRowValues[tableindex + 2] = results[output_index].errors()[index];
                // tableRowValues[curRowValIndex++] = t;
            }
            table_.add_row( tableRowValues );
        }

        if ( print_to_screen_ )
            std::cout << table_ << std::endl;

        if ( save_to_file_ )
        {
            std::ofstream ofs( results_path_ / fs::path( fmt::format( "results-{}.csv", mu.first ) ) );
            VLOG( 2 ) << fmt::format( "[MORTable::update] results file for mu key {}  value {} : {}", mu.first, mu.second, ( results_path_ / fs::path( fmt::format( "results-{}.csv", mu.first ) ) ).string() );
            table_.exportCSV( ofs );
        }
    }
private:    
    MORModels const& models_;
    Table table_;
    bool print_to_screen_;
    bool save_to_file_;
    fs::path results_path_;
};

} // namespace Feel