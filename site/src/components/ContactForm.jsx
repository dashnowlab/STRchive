import { FaPaperPlane } from "react-icons/fa6";
import { Fragment } from "react/jsx-runtime";
import clsx from "clsx";
import { mapValues, startCase, truncate } from "lodash-es";
import { useLocalStorage } from "@reactuses/core";
import { submitContact } from "@/api/contact";
import Form from "@/components/Form";
import TextBox from "@/components/TextBox";
import { repo } from "@/layouts/meta";
import { userAgent } from "@/util/browser";
import { useQuery } from "@/util/hooks";
import { sleep } from "@/util/misc";
import { shortenUrl } from "@/util/string";
import Alert from "./Alert";
import classes from "./Contact.module.css";
import Help from "./Help";
import Link from "./Link";

const ContactForm = () => {
  /** form state, saved to local storage */
  let [name, setName] = useLocalStorage("contact-name", "");
  let [username, setUsername] = useLocalStorage("contact-username", "");
  let [email, setEmail] = useLocalStorage("contact-email", "");
  let [message, setMessage] = useLocalStorage("contact-message", "");

  /** set fallbacks */
  name ||= "";
  username ||= "";
  email ||= "";
  message ||= "";

  /** validate username */
  if (username && username.length > 0)
    username = username.replaceAll(/^@*/g, "@");

  /** extra details to include in report */
  const { pathname, search, hash } = window.location;
  const { browser, engine, os, device, cpu } = userAgent;
  const details = mapValues(
    {
      Page: [pathname + search + hash],
      Browser: [browser.name, browser.version],
      Engine: [engine.name, engine.version],
      OS: [os.name, os.version],
      Device: [device.type, device.model, device.vendor],
      CPU: [cpu.architecture],
    },
    (value) => value.filter(Boolean).join(" "),
  );

  /** message title */
  const title = truncate(
    [name.trim() || username.trim(), message.trim()]
      .filter(Boolean)
      .join(" - "),
    { length: 250 },
  );

  /** message body */
  const body = [{ name, username, email }, details, { message }]
    .map((group) =>
      Object.entries(group)
        .map(([key, value]) => [
          `**${startCase(key)}**`,
          value.trim() ? value.trim() : "\\-",
        ])
        .flat()
        .join("\n"),
    )
    .join("\n\n");

  const {
    query: submit,
    data: response,
    status,
  } = useQuery(async () => {
    try {
      /** test loading spinner */
      await sleep(1000);
      /** will fail because endpoint not set up yet */
      await submitContact(title, body);
    } catch (error) {
      /** prevent error and return fake response */
      return { link: "https://github.com/fake-link" };
    }
  });

  return (
    <Form onSubmit={submit}>
      <div className={clsx("col", classes.form)}>
        <TextBox
          label={
            <>
              Name<Help>Optional. So we know who you are.</Help>
            </>
          }
          placeholder="Your Name"
          value={name}
          onChange={setName}
        />
        <TextBox
          label={
            <>
              GitHub Username
              <Help>
                Optional. So we can tag you in the post and you can follow it.
              </Help>
            </>
          }
          placeholder="@username"
          value={username}
          onChange={setUsername}
        />
        <TextBox
          label={
            <>
              Email
              <Help>Optional. So we can contact you directly if needed.</Help>
            </>
          }
          placeholder="your.name@email.com"
          value={email}
          onChange={setEmail}
        />
        <TextBox
          label="Message"
          placeholder="Comments, questions, issues, etc."
          value={message}
          onChange={setMessage}
          multi
          required
        />
      </div>

      <details>
        <summary>Included details</summary>
        <dl>
          {Object.entries(details).map(([key, value]) => (
            <Fragment key={key}>
              <dt>{key}</dt>
              <dd>{value}</dd>
            </Fragment>
          ))}
        </dl>
      </details>

      <Alert type={status || "info"} style={{ width: 0, minWidth: "100%" }}>
        {status === "" && (
          <>
            This will make a <strong>public</strong> post on{" "}
            <Link to={repo}>our GitHub</Link> with{" "}
            <em>all of the information above</em>. You'll get a link to it once
            it's created.
          </>
        )}
        {startCase(status)}{" "}
        {response && (
          <Link to={response.link}>{shortenUrl(response.link)}</Link>
        )}
      </Alert>

      {status === "" && (
        <button type="submit">
          <FaPaperPlane />
          Submit
        </button>
      )}
    </Form>
  );
};

export default ContactForm;
