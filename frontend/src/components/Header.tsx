import React from 'react';
import { Link, useLocation } from 'react-router-dom';
import { Button } from '@/components/ui/button';
import {
  NavigationMenu,
  NavigationMenuList,
  NavigationMenuItem,
  NavigationMenuLink,
  navigationMenuTriggerStyle,
} from '@/components/ui/navigation-menu';
import './Header.css';
import { DropdownMenu, DropdownMenuContent, DropdownMenuItem, DropdownMenuTrigger } from "@/components/ui/dropdown-menu"

// ✅ CTA Component
const CTA = () => {
  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="default" className="bg-purple-600 hover:bg-purple-700 text-white">
          Get Started
        </Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent className="bg-white shadow-lg rounded-xl p-2">
        <DropdownMenuItem asChild>
          <Link to="/tool1" className="w-full px-3 py-2 rounded-md hover:bg-purple-100 font-semibold">
            SCP
          </Link>
        </DropdownMenuItem>
        <DropdownMenuItem asChild>
          <Link to="/tool2" className="w-full px-3 py-2 rounded-md hover:bg-purple-100 font-semibold">
            TUSM
          </Link>
        </DropdownMenuItem>
      </DropdownMenuContent>
    </DropdownMenu>
  )
}

// ✅ Header Component
const Header = () => {
  const location = useLocation();

  const isActive = (path: string) => location.pathname === path;

  return (
    <header className="biominer-header">
      <div className="container mx-auto px-4 py-4 flex items-center justify-between">
        
        {/* Logo + Name */}
        <div className="flex items-center gap-3">
          <img
            src="/logo1.png" // ✅ better to use /logo1.png directly since it's in public/
            alt="BioMiner Logo"
            className="w-10 h-10 object-contain"
          />
          <div className="biominer-logo">
            <span className="text-3xl font-extrabold text-purple-700">BioMiner</span>
            <span className="text-base text-muted-foreground ml-1">
              Graph Mining for Drug Discovery
            </span>
          </div>
        </div>

        {/* Navigation */}
        <NavigationMenu>
          <NavigationMenuList>
            <NavigationMenuItem>
              <Link to="/">
                <NavigationMenuLink
                  className={`${navigationMenuTriggerStyle()} nav-link ${
                    isActive('/') ? 'active' : ''
                  }`}
                >
                  Home
                </NavigationMenuLink>
              </Link>
            </NavigationMenuItem>

            <NavigationMenuItem>
              <Link to="/algorithms">
                <NavigationMenuLink
                  className={`${navigationMenuTriggerStyle()} nav-link ${
                    isActive('/algorithms') ? 'active' : ''
                  }`}
                >
                  Algorithms
                </NavigationMenuLink>
              </Link>
            </NavigationMenuItem>

            <NavigationMenuItem>
              <Link to="/about">
                <NavigationMenuLink
                  className={`${navigationMenuTriggerStyle()} nav-link ${
                    isActive('/about') ? 'active' : ''
                  }`}
                >
                  About Us
                </NavigationMenuLink>
              </Link>
            </NavigationMenuItem>

            <NavigationMenuItem>
              <Link to="/contact">
                <NavigationMenuLink
                  className={`${navigationMenuTriggerStyle()} nav-link ${
                    isActive('/contact') ? 'active' : ''
                  }`}
                >
                  Contact Us
                </NavigationMenuLink>
              </Link>
            </NavigationMenuItem>
          </NavigationMenuList>
        </NavigationMenu>

        {/* ✅ Dropdown CTA */}
        <CTA />
      </div>
    </header>
  );
};

export default Header;
