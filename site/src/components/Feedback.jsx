import clsx from "clsx";
import { mapValues, startCase, truncate } from "lodash-es";
import { useLocalStorage } from "@reactuses/core";
import { submitFeedback } from "@/api/feedback";
import Form from "@/components/Form";
import TextBox from "@/components/TextBox";
import { userAgent } from "@/util/browser";
import { useQuery } from "@/util/hooks";
import { sleep } from "@/util/misc";
import { shortenUrl } from "@/util/string";
import Alert from "./Alert";
import classes from "./Feedback.module.css";
import Link from "./Link";

const Feedback = () => {
  /** form state, saved to local storage */
  let [name, setName] = useLocalStorage("feedback-name", "");
  let [username, setUsername] = useLocalStorage("feedback-username", "");
  let [email, setEmail] = useLocalStorage("feedback-email", "");
  let [feedback, setFeedback] = useLocalStorage("feedback-body", "");

  /** set fallbacks */
  name ||= "";
  username ||= "";
  email ||= "";
  feedback ||= "";

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

  /** feedback title */
  const title = truncate(
    [name.trim() || username.trim(), feedback.trim()]
      .filter(Boolean)
      .join(" - "),
    { length: 250 },
  );

  /** feedback body */
  const body = [{ name, username, email }, details, { feedback }]
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
      await submitFeedback(title, body);
    } catch (error) {
      /** prevent error and return fake response */
      return { link: "https://github.com/fake-link" };
    }
  });

  return (
    <Form onSubmit={submit}>
      <div className={clsx("col", classes.feedback)}>
        <TextBox label="Name" value={name} onChange={setName} />
        <TextBox label="Username" value={username} onChange={setUsername} />
        <TextBox label="Email" value={email} onChange={setEmail} />
        <TextBox
          label="Feedback"
          value={feedback}
          onChange={setFeedback}
          multi
        />
      </div>

      {status === "" && <button type="submit">Submit</button>}

      {status !== "" && (
        <Alert type={status}>
          {startCase(status)}{" "}
          {response && (
            <Link to={response.link}>{shortenUrl(response.link)}</Link>
          )}
        </Alert>
      )}
    </Form>
  );
};

export default Feedback;
