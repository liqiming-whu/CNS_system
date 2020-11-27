import React from 'react';
import './App.css';
import ArticleList from "./ArticleList";
import {
  Switch,
  Link,
  Route
} from "react-router-dom";
import About from "./About";
import ArticleDetail from "./ArticleDetail";

const Home = () =>
  <ArticleList />

function App() {
  return (
    <div className="App">
      <header>
        <nav>
          <ul>
            <li>
              <Link to="/">Home</Link>
            </li>
            <li>
              <Link to="/about">About</Link>
            </li>
          </ul>
        </nav>
      </header>

      <hr />

      <Switch>
        <Route path="/about">
          <About />
        </Route>
        <Route path="/articles/:articleId">
          <ArticleDetail />
        </Route>
        <Route exact path="/">
          <Home />
        </Route>
      </Switch>
    </div>
  );
}

export default App;
