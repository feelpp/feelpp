find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::shared_ptr/std::shared_ptr/g' file
find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::make_shared/std::make_shared/g' file
find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::enable_shared/std::enable_shared/g' file
find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::weak/std::weak/g' file
find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::dynamic_pointer/std::dynamic_pointer/g' file
find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::static_pointer/std::static_pointer/g' file
find . -name "*.[ch]pp" | xargs -Ifile sed -i.bak 's/boost::const_pointer/std::const_pointer/g' file

