import React, {Component} from "react";
import axios from 'axios';

axios.defaults.withCredentials = true;

class ArticleList extends Component {
    constructor(props) {
        super(props);
        this.state = {
          articleList: []
        }
    }

    componentDidMount() {
      axios.get('/articles/')
      .then(response => {
        this.setState({
          articleList: response.data
        })
      }).catch(err => {
        console.log(err)
      });
    }

    render() {
        return (
            <div className="ArticleList">
            {this.state.articleList.map(item =>
              <div key={item.pmid}>
                <h4>{item.title}</h4>
                <p>
                  <strong>{item.abstract}</strong>
                  <br/>
                  <em>{item.pubdate}</em>
                  <em>{item.topic}</em>
                </p>
              </div>
            )}
          </div>
        );
    }
}

export default ArticleList;