#pragma once


namespace Feel
{
template <typename ToolboxType,typename ... Ts>
int
execute( std::shared_ptr<ToolboxType>& toolbox, Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    int verbose = args.get_else( _verbose, 2 );
    bool save = args.get_else( _save, true );

    toolbox->init();

    if ( verbose > 1 )
        toolbox->printAndSaveInfo();

    if ( toolbox->isStationary() )
    {
        toolbox->solve();
        toolbox->exportResults(/*save*/);
    }
    else
    {
        if ( !toolbox->doRestart() )
            toolbox->exportResults(toolbox->timeInitial()/*, save*/);

        for ( toolbox->startTimeStep() ; !toolbox->timeStepBase()->isFinished(); toolbox->updateTimeStep() )
        {
            if (verbose > 1 && toolbox->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << toolbox->time() << "s \n";
                std::cout << "============================================================\n";
            }

            toolbox->solve();
            toolbox->exportResults(/*save*/);
        }
    }
    return !toolbox->checkResults();
}

} // namespace Feel
