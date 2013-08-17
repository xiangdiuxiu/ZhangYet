Summary:Modified plink 
Name: mplink
Version: 1.0
Release: 1
License: GPL
Group: Development/Tools
URL: https://github.com/
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root

%description

%prep
%setup -q

%build
make
%install
rm -rf $RPM_BUILD_ROOT
RPM_INSTALL_ROOT=$RPM_BUILD_ROOT make install
%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
/usr/local/bin/mplink
%doc


%changelog
* Fri Aug 16 2013 ZhangYet <dante@localhost.localdomain> - 
- Initial build.

